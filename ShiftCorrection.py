#!/usr/local/bin python3

from gi.repository import Gtk
from math import *
import h5py
import hdf5plugin
import numpy as np
from scipy.signal import fftconvolve
from skimage.filters import threshold_otsu
from silx.image import sift
from readMDA import *
from matplotlib.figure import Figure
from matplotlib import colors, patches
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
import itertools


def Display2Data(Axe, x, y):
    """display coordinate to data coordinate"""

    return Axe.transData.inverted().transform(np.array([(x,y)]))[0]


def cross_correlation(im1, im2, mode):

    if mode == "1":
        im_1 = im1-threshold_otsu(im1)
        im_2 = im2-threshold_otsu(im2)
    elif mode == "2":
        im_1 = im1-im1.mean()
        im_2 = im2-im2.mean()
    ccl = fftconvolve(im_1, im_2[::-1,::-1], mode='same')
    dy, dx = np.unravel_index(np.argmax(ccl), im2.shape)
    dy -= im2.shape[0]/2.
    dx -= im2.shape[1]/2.
    return int(round(dy,0)), int(round(dx,0))


def calculate_shift(i_range):

    global xlims, ylims
    
    for i in i_range:
        xlims[i,0] = data[i].sum(0).nonzero()[0].min()
        xlims[i,1] = data[i].sum(0).nonzero()[0].max()
        ylims[i,0] = data[i].sum(1).nonzero()[0].min()
        ylims[i,1] = data[i].sum(1).nonzero()[0].max()
        figs[i].label1.set_text("dY {0} dX {1}       ".format(ylims[i,0]-dim0_max, xlims[i,0]-dim1_max))


def sift_correction(im1, im2, mode):

    if mode == '4':
        im1 = np.log(im1+1)
        im2 = np.log(im2+1)
    sift_ocl = sift.SiftPlan(template=im1)
    mp = sift.MatchPlan()
    match = mp(sift_ocl(im1), sift_ocl(im2))
    dx = np.median(match[:,1].x-match[:,0].x)
    dy = np.median(match[:,1].y-match[:,0].y) 
    return int(round(-dy,0)), int(round(-dx,0)) 


class Figure_with_Projections:

    def on_hover(self, event):
        
        x, y = map(lambda x: int(round(x, 0)), (Display2Data(self.image_ax, event.x, event.y)))

        if Image_Panning:
            datatmp = np.roll(np.copy(data[self.i]), (y-self.y0, x-self.x0), axis=(0,1))
            self.image_im.set_array(datatmp)
            self.canvas.draw()
        elif Image_Zooming:
            if self.non_animated_background != None:
                # restore the clean slate background
                self.canvas.restore_region(self.non_animated_background)
                self.rectangle.set_width(x-self.x0)
                self.rectangle.set_height(y-self.y0)
                self.image_ax.draw_artist(self.rectangle)
                self.canvas.blit(self.image_ax.bbox)
            else:
                self.rectangle =  patches.Rectangle((self.x0, self.y0), width = x-self.x0, height = y - self.y0,\
                                                    alpha = 1, edgecolor = 'k', fill = False, linewidth= 2, animated = True)
                self.image_ax.add_patch(self.rectangle)
                self.canvas.draw()
                self.non_animated_background = self.canvas.copy_from_bbox(self.image_ax.bbox)
 

    def on_click(self, event):

        global Image_Panning, Image_Zooming, lines

        if event.button == 1 or event.button == 3:
            x, y = map(lambda x: int(round(x, 0)), (Display2Data(self.image_ax, event.x, event.y)))
            if Image_Zoomed and event.button == 1:
                try:
                    for line in lines:
                        line.remove()
                except:
                    pass
                lines = []
                for i in range(len(figs)):
                    lines += [figs[i].image_ax.axhline(y, color='k')]
                    lines += [figs[i].image_ax.axvline(x, color='k')]
                    figs[i].canvas.draw()
                    figs[i].label1.set_text("Z {0}       ".format(data[i,y,x]))
                #for f in figs:
                #    lines += [f.image_ax.axhline(y, color='k')]
                #    lines += [f.image_ax.axvline(x, color='k')]
                #    f.canvas.draw()
            else:
                if event.button == 1:
                    Image_Zooming = True
                elif event.button == 3:
                    Image_Panning = True                
                self.x0 = x
                self.y0 = y
                event.canvas.mpl_disconnect(self.on_click_handler)
                self.on_release_handler = self.canvas.mpl_connect('button_release_event', self.on_release)
                

    def on_release(self, event):

        global Image_Panning, Image_Zooming, Image_Zoomed
        
        x, y = map(lambda x: int(round(x, 0)), (Display2Data(self.image_ax, event.x, event.y)))

        event.canvas.mpl_disconnect(self.on_release_handler)
        self.on_click_handler = self.canvas.mpl_connect('button_press_event', self.on_click)

        xmin, xmax = min(x,self.x0), max(x,self.x0)
        ymin, ymax = min(y,self.y0), max(y,self.y0)

        if Image_Zooming:
            try:
                self.rectangle.remove()
            except:
                pass
            self.non_animated_background =  None
            Image_Zooming = False
            Image_Zoomed = True
            for f in figs:
                f.image_ax.set_xlim([xmin-0.5,xmax-0.5])
                f.image_ax.set_ylim([ymin-0.5,ymax-0.5])
                f.canvas.draw()
        elif Image_Panning:
            Image_Panning = False
            data[self.i] = np.roll(data[self.i], (y-self.y0, x-self.x0), axis=(0,1))
            self.image_im.set_array(data[self.i])
            self.canvas.draw()
            calculate_shift([self.i])

            
    def on_press(self, event):

        global xlims, ylims, Image_Zoomed, data, colorscale, theta, cm1
        
        if event.key == "r":
            for f in figs:
                f.image_ax.set_xlim(dim1_max-0.5, dim1_max*2.-0.5)
                f.image_ax.set_ylim(dim0_max-0.5, dim0_max*2.-0.5)
                f.canvas.draw()
            Image_Zoomed = False
        elif event.key == "s":
            calculate_shift(range(data.shape[0]))
            xmin = int(xlims[:,0].max())
            xmax = int(xlims[:,1].min())
            ymin = int(ylims[:,0].max())
            ymax = int(ylims[:,1].min())
            ndim = np.array(list(map(np.shape,dataraw)))
            np.savetxt(shift_file, np.vstack((scannum, theta*1000, ymin-ylims[:,0], ymax-ylims[:,0], xmin-xlims[:,0],xmax-xlims[:,0], ndim[:,0], ndim[:,1], ystep, xstep)).T, fmt='%7d %7d %7d %7d %7d %7d %7d %7d %7d %7d', header=' scan   theta    ylow   yhigh    xlow   xhigh    dimy    dimx   ystep   xstep')
            print ("Shift correction saved!")
        elif event.key == "l":
            try:
                data *= 0
                shifts = np.loadtxt(shift_file, dtype=int)
                for i in range(data.shape[0]):
                    data[i,dim0_max:dim0_max+dataraw[i].shape[0],dim1_max:dim1_max+dataraw[i].shape[1]] = dataraw[i]
                    ymin = int(shifts[:,2].min())
                    xmin = int(shifts[:,4].min())
                    data[i] = np.roll(data[i],(int(ymin-shifts[i,2]), int(xmin-shifts[i,4])),axis=(0,1))
                    figs[i].image_im.set_array(data[i])
                    figs[i].canvas.draw()
                calculate_shift(range(data.shape[0]))
                theta = shifts[:,1]/1000. # also read in theta in case you have overwritten the default values in the saved files
            except Exception as e:
                print(e)
        elif event.key in ["1","2"]:
            xmin,xmax = figs[0].image_ax.get_xlim()
            ymin,ymax = figs[0].image_ax.get_ylim()
            xmin, xmax, ymin, ymax = int(xmin+0.5), int(xmax-0.5), int(ymin+0.5), int(ymax-0.5)
            for i in range(data.shape[0]-1,0,-1):
                data_1 = data[i-1,ymin:ymax+1,xmin:xmax+1]
                data_2 = data[i,ymin:ymax+1,xmin:xmax+1]
                dy, dx = cross_correlation(data_1, data_2, mode=event.key)
                for j in range(i, data.shape[0]):
                    data[j] = np.roll(data[j], (dy, dx), axis=(0,1))
            for i in range(1,data.shape[0]):
                figs[i].image_im.set_array(data[i])
                figs[i].canvas.draw()
            calculate_shift(range(1,data.shape[0]))
        elif event.key in ["3","4"]:
            xmin,xmax = figs[0].image_ax.get_xlim()
            ymin,ymax = figs[0].image_ax.get_ylim()
            xmin, xmax, ymin, ymax = int(xmin+0.5), int(xmax-0.5), int(ymin+0.5), int(ymax-0.5)
            for i in range(data.shape[0]-1,0,-1):
                data_1 = data[i-1,ymin:ymax+1,xmin:xmax+1]
                data_2 = data[i,ymin:ymax+1,xmin:xmax+1]
                dy, dx = sift_correction(data_1, data_2, mode=event.key)
                for j in range(i, data.shape[0]):
                    data[j] = np.roll(data[j], (dy, dx), axis=(0,1))
            for i in range(1,data.shape[0]):
                figs[i].image_im.set_array(data[i])
                figs[i].canvas.draw()
            calculate_shift(range(1,data.shape[0]))
        elif event.key == "u":
            calculate_shift(range(data.shape[0]))
            for f in figs:
                f.image_ax.set_xlim(xlims[:,0].max()-0.5, xlims[:,1].min()+0.5)
                f.image_ax.set_ylim(ylims[:,0].max()-0.5, ylims[:,1].min()+0.5)
                f.canvas.draw()
        elif event.key == "c":
            colorscale = (colorscale+1)%4
            for i in range(data.shape[0]):
                if colorscale == 0:
                    figs[i].image_im.set_norm(colors.Normalize(dataraw_min, dataraw_max))
                elif colorscale == 1:
                    figs[i].image_im.set_norm(colors.Normalize(dataraw[i].min(), dataraw[i].max()))
                elif colorscale == 2:
                    figs[i].image_im.set_norm(colors.LogNorm(max(1,dataraw_min), dataraw_max))
                elif colorscale == 3:
                    figs[i].image_im.set_norm(colors.LogNorm(max(1,dataraw[i].min()), dataraw[i].max()))
                figs[i].canvas.draw()
        elif event.key == "m":
            cm1 = next(cm0)
            for f in figs:
                f.image_im.set_cmap(cm1)
                f.canvas.draw()
        

                
            
    def on_enter(self, event):

        self.canvas.grab_focus()

            
    def __init__(self, index):
        
        self.i = index
        self.non_animated_background = None
        self.fig = Figure()
        self.image_ax = self.fig.add_axes([0,0,1,1])
        self.image_ax.set_axis_off()
        self.image_im = self.image_ax.imshow(data[index], cmap=cm1, norm=colors.Normalize(dataraw_min,dataraw_max), origin="lower", interpolation = "nearest", aspect="equal")
        self.image_ax.set_xlim(dim1_max-0.5, dim1_max*2-0.5)
        self.image_ax.set_ylim(dim0_max-0.5, dim0_max*2-0.5)
        self.canvas = FigureCanvas(self.fig)
        self.on_hover_handler = self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.on_click_handler = self.canvas.mpl_connect('button_press_event', self.on_click)
        self.on_press_handler = self.canvas.mpl_connect('key_press_event', self.on_press)
        self.on_enter_handler = self.canvas.mpl_connect('figure_enter_event', self.on_enter)
        self.vbox = Gtk.VBox()
        self.hbox = Gtk.HBox()
        label0 = Gtk.Label("        Scan {0} Theta {1:.2f}".format(scannum[self.i], theta[self.i]))
        self.label1 = Gtk.Label()
        self.hbox.pack_start(label0, False, False, 0)
        self.hbox.pack_end(self.label1, False, False, 0)
        self.vbox.pack_start(self.hbox, False, False, 0)
        self.vbox.pack_start(self.canvas, True, True, 0)
        

class MainWindow():

    def MainWindow_Destroy(self, widget, dude): 

        dude.ShiftCorrection = None
        self.win.destroy()
    
    def __init__(self, dude):

        global dataraw, data, scannum, theta, dim0_max, dim1_max, xlims, ylims
        global figs, shift_file, Image_Panning, Image_Zooming, Image_Zoomed
        global colorscale, xstep, ystep, dataraw_max, dataraw_min
        global cm0, cm1

        dataraw = []
        dataraw_max = 0
        dataraw_min = 1e50
        scannum = []
        theta = []
        figs = []
        dim0_max = 0
        dim1_max = 0
        Image_Panning = False
        Image_Zooming = False
        Image_Zoomed = False
        colorscale = 0
        xstep = []
        ystep = []
        shift_file= os.path.join(os.path.abspath(os.path.join(dude.MDA_folder, os.pardir)),"Analysis", "shifts.txt")
        mdafiles = []

        cm0 = itertools.cycle(["jet", "coolwarm", "viridis", "gray"])
        cm1 = next(cm0)

        for row in dude.MDA_File_ListStore: 
            if row[14] == True: # row[14] is the multiselect toggle button
                scannum += [int(row[0])]
                mdafiles += [row[13]]

        if len(mdafiles) == 0:
            shifts = np.loadtxt(shift_file, dtype=int)
            scannum = shifts[:,0].astype(int).tolist()
            for scan in scannum:
                mdafiles += [os.path.join(dude.MDA_folder, "26idbSOFT_{0:04d}.mda".format(scan))]

        for i in range(len(mdafiles)):
            mda = readMDA(mdafiles[i], verbose=0) # row[12] contains the path of the mda file
            ystep += [np.int(np.round((mda[1].p[0].data[1]-mda[1].p[0].data[0])*1000,0))]
            xstep += [np.int(np.round((mda[2].p[0].data[0][1]-mda[2].p[0].data[0][0])*1000,0))]
            # below it is assumed that you only use this on 2D data
            datatmp =  np.array(mda[2].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
            if "eiger" in mda[2].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].name and dude.dirty_fix:
                datatmp2 = np.copy(datatmp).flatten()
                datatmp2[:-1] = datatmp2[1:]
                datatmp = datatmp2.reshape(datatmp.shape)
            dim0, dim1 = datatmp.shape
            dataraw += [datatmp]
            dataraw_max = max(dataraw_max, datatmp.max())
            dataraw_min = min(dataraw_min, datatmp.min())
            h5_filename = [f.name for f in os.scandir(dude.h5_folder) if "scan_{0}".format(scannum[i]) in f.name]
            if len(h5_filename):
                h5_filename = os.path.join(dude.h5_folder, h5_filename[0])
                #if os.path.exists(h5_filename):
                with h5py.File(h5_filename, 'r') as h5:
                    theta += [h5["/entry/instrument/26-ID-C/SAMPLE THETA"][()]]
                #h5.close()
            else:
                for ii in range(mda[2].nd):
                    if mda[2].d[ii].name == "atto2:PIC867:1:m1.RBV":
                        break
                theta += [mda[2].d[ii].data[0][0]]
            # this is in case the data acquisition had been aborted

            dim0_max = max(dim0_max, dim0)
            dim1_max = max(dim1_max, dim1)
                 
        #dataraw = [x for _,x in sorted(zip(theta,dataraw))]
        #scannum = [x for _,x in sorted(zip(theta,scannum))]
        #theta.sort()

        # sometimes attoz and attox got swapped, this results in negative theta
        theta = np.array(theta)
        theta[theta<0] += 90

        ndata = len(theta)
        data = np.zeros((ndata, dim0_max*3, dim1_max*3), dtype=np.float64)
        xlims = np.zeros((ndata,2), dtype=np.int32)
        ylims = np.zeros((ndata,2), dtype=np.int32)
        
        nrows = int(sqrt(ndata))
        ncols = int(ndata/nrows) if not ndata%nrows else int(ndata/nrows)+1

        table = Gtk.Table(n_rows = nrows, n_columns = ncols, homogeneous = True)
        for i in range(ndata):
            ylims[i] = dim0_max, dim0_max+dataraw[i].shape[0]-1
            xlims[i] = dim1_max, dim1_max+dataraw[i].shape[1]-1
            data[i,ylims[i,0]:ylims[i,1]+1,xlims[i,0]:xlims[i,1]+1] = dataraw[i]
            figs += [Figure_with_Projections(i)]
            table.attach(figs[-1].vbox, i%ncols,i%ncols+1,i/ncols,i/ncols+1)
        
        instruction_label = Gtk.Label("[left mouse: zoom in / add crosshair] [right mouse: manual shift] [r: reset zoom] [1: xcorr-otsu, 2: xcorr-mean, 3: sift linear, 4: sift log] [u: usable area] [s: save] [l: load] [c: change colorscale] [m: change colormap]")

        table.set_row_spacings(3)  
        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        Main_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Main_VBox.set_border_width(3)
        Main_VBox.pack_start(instruction_label, False, False, 0)
        Main_VBox.pack_start(table, True, True, 0)

        self.win.add(Main_VBox)
        self.win.set_title("Shift Correction")
        self.win.maximize()
        self.win.show_all()
 
