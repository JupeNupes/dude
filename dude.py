#!/usr/local/bin python3
# Standard library and third-party imports for GUI, data processing, and visualization

"""dude.py: Just Another Data Viewer"""

import os
import sys
import time
from readMDA import *
from math import *
import numpy as np
from multiprocessing import Pool, cpu_count
import multiprocessing
multiprocessing.set_start_method('fork', force=True)  # Required on macOS to prevent child processes from re-importing this module
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
from matplotlib import pyplot, colors, patches, font_manager
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import fabio
import h5py
import hdf5plugin
import netCDF4
from scipy import sparse
from misc_dude import *
import Shortcuts
import PowderHelper
import ShowMetadata
import ShiftCorrection
import Analyze5D
try:
    import PtychoLib
except:
    pass
else:
    import ProbeLibrary
    
try:
    import gspread
except:
    pass
else:
    import uuid
    from googleapiclient.discovery import build
    from google_auth_oauthlib.flow import InstalledAppFlow
    from google.auth.transport.requests import Request
    from google.oauth2.credentials import Credentials
    from googleapiclient.http import MediaFileUpload    
    from oauth2client.service_account import ServiceAccountCredentials
    from gspread_formatting import set_row_height

class MyMainWindow:

    #------------------------Add-ons-------------------------------#
    # These methods manage sub-windows for additional functionality (Add-ons)
    # They check if a window already exists and present it, otherwise they initialize it.

    def Addon_PtychoLib(self, widget):
        
        try:
            self.PtychoLib.win.present()
        except:
            self.PtychoLib = PtychoLib.MainWindow(self)

    def Addon_Shortcuts(self, widget):
        
        try:
            self.Shortcuts.win.present()
        except:
            self.Shortcuts = Shortcuts.MainWindow(self)

    def Addon_ProbeLibrary(self, widget):
        
        try:
            self.ProbeLibrary.win.present()
        except:
            self.ProbeLibrary = ProbeLibrary.MainWindow(self)

    def Addon_PowderHelper(self, widget):
        
        try:
            self.PowderHelper.win.present()
        except:
            self.PowderHelper = PowderHelper.MainWindow(self)

    def Addon_ShiftCorrection(self, widget):
        
        try:
            self.ShiftCorrection.win.present()
        except:
            self.ShiftCorrection = ShiftCorrection.MainWindow(self)

    def Addon_Analyze5D(self, widget):
        
        try:
            self.Analyze5D.win.present()
        except:
            self.Analyze5D = Analyze5D.MainWindow(self)

    def Addon_ShowMetadata(self, widget):
        
        try:
            self.ShowMetadata.win.present()
        except:
            self.ShowMetadata = ShowMetadata.MainWindow(self)

    #-------------------------miscellaneous------------------------#
    # Utility methods for console commands and canvas drawing


    def Console_Command_Sent(self, widget):
        # Handles commands entered in the GUI console (clear, history, or arbitrary python code)


        if self.Console_Entry.get_text() == 'clear':
            widget.get_completion().get_model().clear()
        elif self.Console_Entry.get_text() == 'history':
            for command in widget.get_completion().get_model():
                print(command[0])
        else:
            try:
                exec(self.Console_Entry.get_text())
            except Exception as e:
                print(e)
            else:
                widget.get_completion().get_model().append([self.Console_Entry.get_text()])
                self.console_cmd_index = len(widget.get_completion().get_model())
                self.Console_Entry.set_text('')


    def Image_Canvas_Draw_Rectangle(self, xmin, xmax, ymin, ymax, color, animated = False):
        # Helper to draw a rectangle on the detector image axe (used for ROIs and Zooming)


        Rectangle =  patches.Rectangle((xmin, ymin), width = xmax-xmin, height = ymax - ymin,\
        alpha = 1, edgecolor = color, fill = False, linewidth= 3, animated = animated)
        self.Image_Axe.add_patch(Rectangle)
        return Rectangle

    #------------detector image related callbacks--------------#
    # Methods for loading and interacting with detector images


    def LoadAllImages(self, widget):
        # Loads all detector images for the current scan into memory.
        # Supports HDF5 (Eiger) and TIFF stacks, with options for sparse matrix storage and masking.

        
        t0 = time.time()
        if self.eiger_enabled:
            dd = self.h5["/entry/instrument/detector/data"]
            nd0 = dd.maxshape[0]-dd.shape[0] # maxshape is the intended shape, before interruption of scans
            if self.pump_probe:
                nd0 = int(nd0/2)
            if self.sparse_enabled:
                self.image = None
                nd1 = dd.shape[0]
                nd3 = 0
                nblock = 1000
                for i in range(int(nd1/nblock)+(nd1%nblock>0)):
                    nd2 = nblock if nd1>nblock else nd1
                    nd1 -= nblock
                    if not self.pump_probe:
                        ddd = dd[nd3:nd3+nd2][()].reshape(nd2, 1062*1028)
                        ddd[ddd>2**31]=0
                    else:
                        ddd = np.zeros((int(nd2/2), 1062, 1028))
                        tmp = dd[nd3:nd3+nd2,-512:,:][()].reshape(int(nd2/2), 2, 512, 1028)
                        ddd[:,:512,:] = tmp[:,0]
                        ddd[:,-512:,:] = tmp[:,1]
                        ddd = ddd.reshape(int(nd2/2), 1062*1028)
                        ddd[ddd>65530]=0
                    if i == 0:
                        self.image = sparse.csr_matrix(ddd)
                    else:
                        self.image = sparse.vstack((self.image, sparse.csr_matrix(ddd)))
                    nd3 += nblock
                if nd0:
                    self.image = sparse.vstack((self.image, sparse.csr_matrix(np.zeros((nd0,1062*1028)))))
                #self.image = None
                #tmpdata = self.h5["/entry/instrument/detector/data"][()]
                #tmpdata[tmpdata>2**13]=0
                #if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                #    tmpdata *= self.custom_mask
                #self.image = sparse.csr_matrix(tmpdata.reshape(self.h5["/entry/instrument/detector/data"].shape[0], 1062*1028))
                #tmpdata = None
            else:
                self.image = dd[()]
                self.image[self.image>2**31] = 0
                if nd0:
                    self.image = np.vstack((self.image, np.zeros((nd0, 1062, 1028))))
                if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                    self.image *= self.custom_mask
        else:
            image_list = list(map(lambda x:os.path.join(self.image_path,x), self.image_list))
            self.image = np.zeros((self.ntiff,self.dimY, self.dimX), dtype="int32") # edist 20211031 from int16
            pool = Pool(processes=cpu_count()) 
            for i, tmparr in enumerate(pool.imap(image_loader, image_list)):
                self.image[i] = tmparr
            pool.close()
            self.image[self.image<0] = 0
            if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                self.image *= self.custom_mask
        print (time.time()-t0)
        self.AllImagesLoaded = True
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(True)


    def Image_Canvas_Button_Released(self, event):
        # Handles the completion of mouse interactions (Zooming, Panning, ROI drawing) on the detector image.


        self.Image_xend, self.Image_yend = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))
       
        self.Image_Canvas.mpl_disconnect(self.Image_Canvas_Button_Release_Event)
        self.Image_Canvas_Button_Press_Event = self.Image_Canvas.mpl_connect('button_press_event', self.Image_Canvas_Button_Pressed)

        if self.Image_Rectangle_Drawing:
            self.Image_Rectangle_Drawing = False
            self.non_animated_background = None
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            xmin = int(min(self.Image_xend, self.Image_xstart))
            xmax = int(max(self.Image_xend, self.Image_xstart))
            ymin = int(min(self.Image_yend, self.Image_ystart))
            ymax = int(max(self.Image_yend, self.Image_ystart))
            
            self.Image_ShowROI_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b')
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(xmin)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(xmax)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(ymin)
            self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
            self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(ymax)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
        elif self.Image_Zooming:
            self.Image_Zooming = False
            self.Image_Zoomed = True
            self.Image_ZoomOut_Button.set_sensitive(True)
            self.non_animated_background = None
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            xmin = int(min(self.Image_xend, self.Image_xstart))
            xmax = int(max(self.Image_xend, self.Image_xstart))
            ymin = int(min(self.Image_yend, self.Image_ystart))
            ymax = int(max(self.Image_yend, self.Image_ystart))
            dim_y, dim_x = self.Image_Image.get_array().shape
            dx = xmax - xmin
            dy = ymax - ymin
            if 1.*dim_y/dim_x > dy/dx:
                self.Image_Axe.set_xlim(xmin-0.5, xmax-0.5)
                self.Image_Axe.set_ylim((ymin+ymax)/2.+.5*dx/dim_x*dim_y-0.5, (ymin+ymax)/2.-.5*dx/dim_x*dim_y-0.5)
            else:
                self.Image_Axe.set_ylim(ymax-0.5, ymin-0.5)
                self.Image_Axe.set_xlim((xmin+xmax)/2.-.5*dy/dim_y*dim_x-0.5, (xmin+xmax)/2.+.5*dy/dim_y*dim_x-0.5)
            self.Image_Canvas.draw()
        elif self.Image_Panning:
            self.Image_Panning = False
            self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.TCROSS))
       

    def Image_Canvas_Button_Scrolled(self, event):
        # Implements zooming in/out on the detector image using the mouse scroll wheel.


        xdata, ydata = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))
        old_x_range = self.Image_Axe.get_xlim()[1]-self.Image_Axe.get_xlim()[0]
        old_y_range = self.Image_Axe.get_ylim()[1]-self.Image_Axe.get_ylim()[0]
        range_x_new = old_x_range*(1-0.1*event.step)
        range_y_new = old_y_range*(1-0.1*event.step)

        self.Image_Axe.set_xlim(xdata-range_x_new*(xdata-self.Image_Axe.get_xlim()[0])/old_x_range, xdata+range_x_new*(self.Image_Axe.get_xlim()[1]-xdata)/old_x_range)
        self.Image_Axe.set_ylim(ydata-range_y_new*(ydata-self.Image_Axe.get_ylim()[0])/old_y_range, ydata+range_y_new*(self.Image_Axe.get_ylim()[1]-ydata)/old_y_range)
        self.Image_Canvas.draw()
        self.Image_ZoomOut_Button.set_sensitive(True)
        self.Image_Zoomed = True
                

    def Image_Canvas_Figure_Entered(self, event):

        #Otherwise button_press_event.key produces none
        self.Image_Canvas.grab_focus()


    def Image_Canvas_Mouse_Hover(self, event):
        # Updates status labels and interactive overlays as the mouse moves over the detector image.


        self.Image_xend, self.Image_yend = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Image_Axe, event.x, event.y))))
        try:
            dim_y, dim_x = self.Image_Image.get_array().shape
        except:
            return

        if (self.Image_yend < dim_y-0.5) and (self.Image_yend > -.5) and (self.Image_xend > -0.5) and (self.Image_xend < dim_x-0.5):
            if self.show_angle:
                self.Image_Value_Label.set_text('[2The={0:.2f}]  [Gam={1:.2f}]  [Z={2}]'.format(self.twotheta[self.Image_yend, self.Image_xend], self.gamma[self.Image_yend, self.Image_xend], int(self.Image_Image.get_array().data[self.Image_yend, self.Image_xend])))
                # since i don't know which version of matplotlib, somehow get_array() returns a masked array hence the new .data
            else:
                self.Image_Value_Label.set_text('[X = {0:d}]  [Y = {1:d}]  [Z = {2}]'.format(self.Image_xend, self.Image_yend, int(self.Image_Image.get_array().data[self.Image_yend, self.Image_xend])))

            if self.Image_Rectangle_Drawing or self.Image_Zooming:
                if self.non_animated_background != None:
                    # restore the clean slate background
                    self.Image_Canvas.restore_region(self.non_animated_background)
                    if self.Image_xstart > self.Image_xend: #modify the starting point
                        self.Image_ShowROI_Rectangle.set_x(self.Image_xend)
                    self.Image_ShowROI_Rectangle.set_width(abs(self.Image_xend-self.Image_xstart))
                    if self.Image_ystart > self.Image_yend: #modify the starting point
                        self.Image_ShowROI_Rectangle.set_y(self.Image_yend)
                    self.Image_ShowROI_Rectangle.set_height(abs(self.Image_yend-self.Image_ystart))
                    self.Image_Axe.draw_artist(self.Image_ShowROI_Rectangle)
                    self.Image_Canvas.blit(self.Image_Axe.bbox)
                else:
                    xmax = int(max(self.Image_xend, self.Image_xstart))
                    xmin = int(min(self.Image_xend, self.Image_xstart))
                    ymax = int(max(self.Image_yend, self.Image_ystart))
                    ymin = int(min(self.Image_yend, self.Image_ystart))
                    self.Image_ShowROI_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b', animated = True)
                    self.Image_Canvas.draw()
                    self.non_animated_background = self.Image_Canvas.copy_from_bbox(self.Image_Axe.bbox)
            elif self.Image_Panning:
                self.Image_Axe.set_xlim(self.Image_Axe.get_xlim()[0]+self.Image_xstart-self.Image_xend, \
                                      self.Image_Axe.get_xlim()[1]+self.Image_xstart-self.Image_xend)
                self.Image_Axe.set_ylim(self.Image_Axe.get_ylim()[0]+self.Image_ystart-self.Image_yend, \
                                      self.Image_Axe.get_ylim()[1]+self.Image_ystart-self.Image_yend)
                self.Image_xstart, self.Image_ystart = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Image_Axe, event.x, event.y))))
                self.Image_Canvas.draw()
            
    def Image_ZoomOut(self, widget):
        # Resets the detector image view to show the full detector area.


        dim_y, dim_x = self.Image_Image.get_array().shape
        self.Image_Axe.set_xlim(-0.5, dim_x-0.5)
        self.Image_Axe.set_ylim(dim_y-0.5, -0.5)
        self.Image_Canvas.draw()
        self.Image_Zoomed = False
        self.Image_ZoomOut_Button.set_sensitive(False)

    def Image_Canvas_Button_Pressed(self, event):
        # Initiates mouse interactions (Zooming, Panning, ROI drawing) on the detector image.

        
        self.Image_xstart, self.Image_ystart = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))

        if event.button == 2:
            if self.Image_Zoomed:
                dim_y, dim_x = self.Image_Image.get_array().shape
                self.Image_Axe.set_xlim(-0.5, dim_x-0.5)
                self.Image_Axe.set_ylim(dim_y-0.5, -0.5)
                self.Image_Canvas.draw()
                self.Image_Zoomed = False
        else:
            if event.button == 1:
                if self.Image_Zoomed == False:
                    self.Image_Zooming = True
                    try:
                        self.Image_ShowROI_Rectangle.remove()
                    except (AttributeError, ValueError):
                        pass
                else:
                    self.Image_Panning = True
                    self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.FLEUR))
            elif event.button == 3:
                self.Image_Rectangle_Drawing = True
                try:
                    self.Image_ShowROI_Rectangle.remove()
                except (AttributeError, ValueError):
                    pass
            self.Image_Canvas.mpl_disconnect(self.Image_Canvas_Button_Press_Event)
            self.Image_Canvas_Button_Release_Event = self.Image_Canvas.mpl_connect('button_release_event', self.Image_Canvas_Button_Released)


    def Image_AutoScale_toggled(self, widget):
        # Toggles between automatic and manual intensity scaling for the detector image display.


        if widget.get_active():
            data = self.Image_Image.get_array()
            vmin = data.min()
            vmax = data.max()
            self.Image_Vmin_HScale.set_sensitive(False)
            self.Image_Vmax_HScale.set_sensitive(False)
            self.Image_Vmin_HScale_Adjustment.handler_block(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_block(self.Image_Vmax_Changed_Handler)
            self.Image_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Image_Vmin_HScale_Adjustment.set_value(vmin)
            self.Image_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Image_Vmax_HScale_Adjustment.set_value(vmax)
            self.Image_Vmin_HScale_Adjustment.handler_unblock(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_unblock(self.Image_Vmax_Changed_Handler)
            if self.Image_Log_ToggleButton.get_active():
                self.Image_Image.set_norm(colors.LogNorm(vmin = vmin+.1, vmax = vmax+.1))
            else:
                self.Image_Image.set_norm(colors.Normalize(vmin = vmin, vmax = vmax))
            self.Image_Canvas.draw()
        else:
            self.Image_Vmin_HScale.set_sensitive(True)
            self.Image_Vmax_HScale.set_sensitive(True)


    def Image_Vscale_Changed(self, widget):
        # Updates the intensity normalization (linear or log) based on manual slider values.

  
        self.Image_Vmax_HScale_Adjustment.set_lower(self.Image_Vmin_HScale_Adjustment.get_value()+1)
        self.Image_Vmin_HScale_Adjustment.set_upper(self.Image_Vmax_HScale_Adjustment.get_value()-1)
        if self.Image_Log_ToggleButton.get_active():
            self.Image_Image.norm = colors.LogNorm(vmin = self.Image_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Image_Vmax_HScale_Adjustment.get_value()+.1)
        else:
            self.Image_Image.set_norm(colors.Normalize(vmin = self.Image_Vmin_HScale_Adjustment.get_value(), vmax = self.Image_Vmax_HScale_Adjustment.get_value()))
        self.Image_Canvas.draw()


    def Image_Plot_HScale_Changed(self, widget):
        # Updates the displayed detector image when the scan index slider is moved.


        index = int((widget.get_value()-self.image_index_min)/self.nbin)
        if self.eiger_enabled:
            if self.image_index_max > 0:
                if self.pump_probe:
                    image_data = np.zeros((1062,1028), dtype=np.uint16)
                    image_data[:512] = self.h5["/entry/instrument/detector/data"][int(index*2),-512:]
                    image_data[-512:] = self.h5["/entry/instrument/detector/data"][int(index*2+1), -512:]
                else:
                    image_data = self.h5["/entry/instrument/detector/data"][index].astype(np.int32)
            else:
                return # h5 not available
        else:
            image_data = fabio.open(os.path.join(self.image_path,self.image_list[index])).data
        image_data[image_data<0]=0
        if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
            image_data *= self.custom_mask
        try:
            self.Image_Image.set_array(image_data)
        except:
            self.Image_Image = self.Image_Axe.imshow(image_data, aspect = 'equal', interpolation = 'nearest', cmap = self.cm)
        else:
            pass

        vmin = image_data.min()
        vmax = image_data.max()
        self.Image_Vmin_HScale_Adjustment.set_lower(vmin)
        self.Image_Vmax_HScale_Adjustment.set_upper(vmax)
        if self.Image_AutoScale_ToggleButton.get_active():
            self.Image_Vmin_HScale_Adjustment.handler_block(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_block(self.Image_Vmax_Changed_Handler)
            self.Image_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Image_Vmin_HScale_Adjustment.set_value(vmin)
            self.Image_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Image_Vmax_HScale_Adjustment.set_value(vmax)
            self.Image_Vmin_HScale_Adjustment.handler_unblock(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_unblock(self.Image_Vmax_Changed_Handler)
        if self.Image_Log_ToggleButton.get_active():
            self.Image_Image.set_norm(colors.LogNorm(vmin = self.Image_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Image_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Image_Image.set_norm(colors.Normalize(vmin = self.Image_Vmin_HScale_Adjustment.get_value(), vmax = self.Image_Vmax_HScale_Adjustment.get_value()))
        self.Image_Canvas.draw()

    def XRF_RoI(self, widget, flag):
        # Calculates and plots 1D/2D data for a selected XRF Region of Interest.


        if flag == 1:
            xmin = int(self.XRF1_Vmin_HScale_Adjustment.get_value())
            xmax = int(self.XRF1_Vmax_HScale_Adjustment.get_value())
            if self.xrf1_spec == 0:
                ydata = self.xrf_data[:-1,:,xmin:xmax+1].sum(0).sum(1)
            elif self.xrf1_spec == 8:
                ydata = self.xrf_data[-1,:,xmin:xmax+1].sum(1)
            else:
                ydata = self.xrf_data[self.xrf1_spec-1,:,xmin:xmax+1].sum(1)
        else:
            xmin = int(self.XRF2_Vmin_HScale_Adjustment.get_value())
            xmax = int(self.XRF2_Vmax_HScale_Adjustment.get_value())
            if self.xrf2_spec == 0:
                ydata = self.xrf_data[:-1,:,xmin:xmax+1].sum(0).sum(1)
            elif self.xrf2_spec == 8:
                ydata = self.xrf_data[-1,:,xmin:xmax+1].sum(1)
            else:
                ydata = self.xrf_data[self.xrf2_spec-1,:,xmin:xmax+1].sum(1)
        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        if ndim == 1:
            xdata = np.array(self.data[1].p[0].data)
            self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata[(xdata!=0)*(ydata!=0)])
            self.Plot1D_Canvas.draw()
        elif ndim == 2:
            xdata2 = np.zeros((self.data[0]["dimensions"][0],self.data[0]["dimensions"][1]))
            xdata1 = np.copy(xdata2)+1
            datatmp = np.array(self.data[2].p[0].data)
            xdata2[:datatmp.shape[0]] = datatmp
            xdata1 *= np.array(self.data[1].p[0].data)[:,np.newaxis]
            ydata2 = np.zeros(xdata2.shape).flatten()
            ydata2[:ydata.shape[0]]=ydata
            ydata2 = ydata2.reshape(xdata2.shape)
            self.Scan_Plot2D(xdata1, xdata2, ydata2)
            self.Plot2D_Image.set_norm(colors.Normalize(ydata.min(), ydata.max()))
            self.Plot2D_Canvas.draw()
            

    def XRF_Sum(self, widget):
        # Computes the sum of XRF spectra across all points and updates the XRF display axes.


        self.XRF1_Axe.cla()
        index = int(self.XRF_Plot_HScale_Adjustment.get_value())
        if self.xrf1_spec == 0:
            self.xrf1_data = self.xrf_data[:-1].sum(0).sum(0)
        elif self.xrf1_spec == 8:
            self.xrf1_data = self.xrf_data[-1].sum(0)
        else:
            self.xrf1_data = self.xrf_data[self.xrf1_spec-1].sum(0)
        self.XRF1_Axe.plot(self.xrf1_data)
        xmin = int(self.XRF1_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF1_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF1_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf1_data[xmin:xmax+1].min()
            ymax = self.xrf1_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF1_Axe.set_ylim(ymin-dy, ymax+dy)
        self.XRF2_Axe.cla()
        if self.xrf2_spec == 0:
            self.xrf2_data = self.xrf_data[:-1].sum(0).sum(0)
        elif self.xrf2_spec == 8:
            self.xrf2_data = self.xrf_data[-1].sum(0)
        else:
            self.xrf2_data = self.xrf_data[self.xrf2_spec-1].sum(0)
        self.XRF2_Axe.plot(self.xrf2_data)
        xmin = int(self.XRF2_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF2_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF2_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf2_data[xmin:xmax+1].min()
            ymax = self.xrf2_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF2_Axe.set_ylim(ymin-dy, ymax+dy)
        self.XRF_Canvas.draw()


    def XRF1_Vscale_Changed(self, widget):

        xmin = int(self.XRF1_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF1_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF1_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf1_data[xmin:xmax+1].min()
            ymax = self.xrf1_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF1_Axe.set_ylim(ymin-dy, ymax+dy)
            self.XRF_Canvas.draw()


    def XRF2_Vscale_Changed(self, widget):

        xmin = int(self.XRF2_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF2_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF2_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf2_data[xmin:xmax+1].min()
            ymax = self.xrf2_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF2_Axe.set_ylim(ymin-dy, ymax+dy)
        self.XRF1_Vmin_HScale_Adjustment.handler_block(self.XRF1_Vmin_Changed_Handler)
        self.XRF1_Vmin_HScale_Adjustment.set_value(xmin)
        self.XRF1_Vmin_HScale_Adjustment.handler_unblock(self.XRF1_Vmin_Changed_Handler)
        self.XRF1_Vmax_HScale_Adjustment.handler_block(self.XRF1_Vmax_Changed_Handler)
        self.XRF1_Vmax_HScale_Adjustment.set_value(xmax)
        self.XRF1_Vmax_HScale_Adjustment.handler_unblock(self.XRF1_Vmax_Changed_Handler)
        self.XRF1_Vscale_Changed(None)
       

    def XRF1_Spec_Changed(self, widget, flag):

        self.xrf1_spec = int(flag)
        self.XRF1_Plot()
        self.XRF_Canvas.draw()

    def XRF2_Spec_Changed(self, widget, flag):

        self.xrf2_spec = int(flag)
        self.XRF2_Plot()
        self.XRF_Canvas.draw()

    def XRF_Plot_HScale_Changed(self, widget):

        self.XRF1_Plot()
        self.XRF2_Plot()
        self.XRF_Canvas.draw()

    def XRF1_Plot(self, index=None):

        self.XRF1_Axe.cla()
        if index == None:
            index = int(self.XRF_Plot_HScale_Adjustment.get_value())
        if self.xrf1_spec == 0:
            self.xrf1_data = self.xrf_data[:-1,index].mean(0)
        elif self.xrf1_spec == 8:
            self.xrf1_data = self.xrf_data[-1,index]
        else:
            self.xrf1_data = self.xrf_data[self.xrf1_spec-1,index]
        self.XRF1_Axe.plot(self.xrf1_data)
        xmin = int(self.XRF1_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF1_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF1_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf1_data[xmin:xmax+1].min()
            ymax = self.xrf1_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF1_Axe.set_ylim(ymin-dy, ymax+dy)

    def XRF2_Plot(self, index=None):

        self.XRF2_Axe.cla()
        if index == None:
            index = int(self.XRF_Plot_HScale_Adjustment.get_value())
        if self.xrf2_spec == 0:
            self.xrf2_data = self.xrf_data[:-1,index].mean(0)
        elif self.xrf2_spec == 8:
            self.xrf2_data = self.xrf_data[-1,index]
        else:
            self.xrf2_data = self.xrf_data[self.xrf2_spec-1,index]
        self.XRF2_Axe.plot(self.xrf2_data)
        xmin = int(self.XRF2_Vmin_HScale_Adjustment.get_value())
        xmax = int(self.XRF2_Vmax_HScale_Adjustment.get_value())
        if xmin < xmax:
            self.XRF2_Axe.set_xlim(xmin-1, xmax+2)
            ymin = self.xrf2_data[xmin:xmax+1].min()
            ymax = self.xrf2_data[xmin:xmax+1].max()
            dy = (ymax-ymin)*0.05
            self.XRF2_Axe.set_ylim(ymin-dy, ymax+dy)

    #----------------1D Plot related callbacks----------------#
    # Methods for interacting with 1D scan plots


    def Plot1D_Canvas_Mouse_Hover(self, event):

        xpointer, ypointer = Display2Data(self.Plot1D_Axe, event.x, event.y)
        try:
            index = np.fabs(self.Plot1D_xdata - xpointer).argmin()
        except:
            return
        if (self.xrf_mode and index != self.XRF_Plot_HScale_Adjustment.get_value()) or \
           ((not self.xrf_mode) and (self.image_index_min + index * self.nbin != self.Image_Plot_HScale_Adjustment.get_value())):
            if event.button == 3: 
                try:
                    for hl in self.Plot1D_Hightlight:
                        hl.remove()
                except (AttributeError, ValueError):
                    pass
                x, y = self.Plot1D_xdata[index], self.Plot1D_ydata[index]
                self.Plot1D_Hightlight = self.Plot1D_Axe.plot(x, y, "ro", markersize=10)
                self.Plot1D_Canvas.draw()
                self.Plot1D_Value_Label.set_text('[X = {0:.3f}] [Y = {1:.3f}]'.format(x, y))
                if self.xrf_mode:
                    self.XRF_Plot_HScale_Adjustment.set_value(index)
                else:
                    self.Image_Plot_HScale_Adjustment.set_value(self.image_index_min + index * self.nbin)


    def Plot1D_Vscale_Changed(self, widget):

        self.Plot1D_Axe.set_xscale('log' if self.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        self.Plot1D_Axe.set_yscale('log' if self.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        self.Plot1D_Canvas.draw()


    def Scan_Plot1D(self, xdata, ydata, hold=False):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active() and not hold:
            self.Plot1D_Axe.cla()
        self.Plot1D_Axe.plot(xdata, ydata, "-o")
        self.Plot1D_Axe.set_xscale('log' if self.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        self.Plot1D_Axe.set_yscale('log' if self.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        self.Plot1D_xdata = xdata
        self.Plot1D_ydata = ydata

        
    #----------------2D Plot related callbacks----------------#
    # Methods for interacting with 2D scan plots (maps)

    def Plot2D_Canvas_Mouse_Pressed(self, event):

        global _xstart, _ystart
        _xstart, _ystart = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        self.Plot2D_Canvas.mpl_disconnect(self.Plot2D_Canvas_Mouse_Pressed)
        self.Plot2D_Canvas_Button_Release_Event = self.Plot2D_Canvas.mpl_connect('button_release_event', self.Plot2D_Canvas_Mouse_Released)

    def Plot2D_Canvas_Mouse_Released(self, event):

        _xend, _yend = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        self.ttlist += [[(_xstart+event.xdata)/2, (_ystart+event.ydata)/2]]
        self.Plot2D_Canvas.mpl_disconnect(self.Plot2D_Canvas_Mouse_Released)
        self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Pressed)
        self.Plot2D_Axe.plot(self.ttlist[-1][0], self.ttlist[-1][1], 'r.')
        self.Plot2D_Canvas.draw()
    
    def Plot2D_Canvas_Mouse_Hover(self, event):

        x1, x2 = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        try:
            dim_y, dim_x = self.Plot2D_Image.get_array().shape
        except:
            return
        if (x2 < dim_y-0.5) and (x2 > -.5) and (x1 > -0.5) and (x1 < dim_x-0.5):
            if self.xrf_mode:
                index = int(x1+x2*self.Plot2D_ydata.shape[1])
            else:
                index = int(self.image_index_min+(x1+x2*self.Plot2D_ydata.shape[1])*self.nbin)
            if (self.xrf_mode and index != self.XRF_Plot_HScale_Adjustment.get_value()) or \
               ((not self.xrf_mode) and index != self.Image_Plot_HScale_Adjustment.get_value()):
                if event.button == 3:
                    self.Plot2D_P0_Label.set_text('P0 [X = {0:.3f}] [Y = {1:.3f}] [Z = {2:.3f}]'.format(self.Plot2D_xdata2[x2,x1], self.Plot2D_xdata1[x2,x1], self.Plot2D_ydata[x2,x1]))
                    if self.xrf_mode:
                        self.XRF_Plot_HScale_Adjustment.set_value(index)
                    else:
                        self.Image_Plot_HScale_Adjustment.set_value(index)
                elif event.button == 1:
                    lelements = self.Plot2D_P0_Label.get_text().split()
                    x0, y0 = float(lelements[3][:-1]), float(lelements[6][:-1])
                    self.Plot2D_P1_Label.set_text('P1 [X = {0:.3f}] [Y = {1:.3f}] [dX = {2:.3f}] [dY = {3:.3f}]'.format(self.Plot2D_xdata2[x2,x1], self.Plot2D_xdata1[x2,x1], self.Plot2D_xdata2[x2,x1]-x0, self.Plot2D_xdata1[x2,x1]-y0))
                    if self.xrf_mode:
                        self.XRF1_Plot(index)
                        self.XRF_Canvas.draw()
                #elif event.xdata != None and event.ydata != None:
                #    self.Plot2D_Value_Label.set_text('[X = {0:.3f}]  [Y = {1:.3f}]'.format(event.xdata, event.ydata))
    

    def Plot2D_AutoScale_toggled(self, widget):

        if widget.get_active():
            data = self.Plot2D_Image.get_array()
            vmin = data.min()
            vmax = data.max()
            self.Plot2D_Vmin_HScale.set_sensitive(False)
            self.Plot2D_Vmax_HScale.set_sensitive(False)
            self.Plot2D_Vmin_HScale_Adjustment.handler_block(self.Plot2D_Vmin_Changed_Handler)
            self.Plot2D_Vmax_HScale_Adjustment.handler_block(self.Plot2D_Vmax_Changed_Handler)
            self.Plot2D_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Plot2D_Vmin_HScale_Adjustment.set_value(vmin)
            self.Plot2D_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Plot2D_Vmax_HScale_Adjustment.set_value(vmax)
            self.Plot2D_Vmin_HScale_Adjustment.handler_unblock(self.Plot2D_Vmin_Changed_Handler)
            self.Plot2D_Vmax_HScale_Adjustment.handler_unblock(self.Plot2D_Vmax_Changed_Handler)
            if self.Plot2D_Log_ToggleButton.get_active():
                self.Plot2D_Image.set_norm(colors.LogNorm(vmin = vmin+.1, vmax = vmax+.1))
            else:
                self.Plot2D_Image.set_norm(colors.Normalize(vmin = vmin, vmax = vmax))
            self.Plot2D_Canvas.draw()
        else:
            self.Plot2D_Vmin_HScale.set_sensitive(True)
            self.Plot2D_Vmax_HScale.set_sensitive(True)

    def Plot2D_ScaleBar_toggled(self, widget):
        
        try:
            self.Plot2D_scalebar.remove()
        except:
            pass
        if widget.get_active() and self.data[1].curr_pt > 2:
            um_per_pixel = np.abs(self.data[1].p[0].data[1]-self.data[1].p[0].data[0])
            um_choice = np.array([0.1,0.2,0.5,1,2,5,10,20,30,40,50,100,200,500])
            x_um = um_choice[np.argmin(np.abs(um_choice / (um_per_pixel * self.data[2].npts) - 1/4))]
            x_pixel = int(x_um / um_per_pixel)
            self.Plot2D_scalebar = AnchoredSizeBar(self.Plot2D_Axe.transData,\
                                                   x_pixel, u'{0} μm'.format(x_um), 'lower right', \
                                                   pad=1, label_top = True,\
                                                   color='white',\
                                                   frameon=False,\
                                                   size_vertical=self.Plot2D_xdata2.shape[0]/80., \
                                                   fontproperties = font_manager.FontProperties(size=20))
            self.Plot2D_Axe.add_artist(self.Plot2D_scalebar)
        else:
            pass
        self.Plot2D_Canvas.draw()


    def Plot2D_Vscale_Changed(self, widget):
  
        self.Plot2D_Vmax_HScale_Adjustment.set_lower(self.Plot2D_Vmin_HScale_Adjustment.get_value()+1)
        self.Plot2D_Vmin_HScale_Adjustment.set_upper(self.Plot2D_Vmax_HScale_Adjustment.get_value()-1)
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_Canvas.draw()


    def Scan_PlotSpiral(self, xdata1, xdata2, ydata):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active():
            self.Plot2D_Axe.cla()
        self.Plot2D_Axe.set_axis_off()
        self.Plot2D_Image = self.Plot2D_Axe.scatter(xdata2, xdata1, ydata, s=2, edgecolor='none', cmap=self.cm)
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_xdata1 = xdata1
        self.Plot2D_xdata2 = xdata2
        self.Plot2D_ydata = ydata


    def Scan_Plot2D(self, xdata1, xdata2, ydata):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active():
            self.Plot2D_Axe.cla()
        self.Plot2D_Axe.set_axis_off()
        self.Plot2D_Image = self.Plot2D_Axe.imshow(ydata,  aspect = 'equal', interpolation = 'nearest', cmap=self.cm, origin="lower")
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = max(0.1,self.Plot2D_Vmin_HScale_Adjustment.get_value()), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_xdata1 = xdata1
        self.Plot2D_xdata2 = xdata2
        self.Plot2D_ydata = ydata
        if self.Plot2D_ScaleBar_ToggleButton.get_active():
            self.Plot2D_ScaleBar_toggled(self.Plot2D_ScaleBar_ToggleButton)

    #-------------------main menu callbacks--------------------#
    # Methods triggered by main menu items and global settings


    def FileDialog_Construction(self, widget, flag = 0):
        # Generic file/folder dialog constructor for MDA, Images, or Save/Open actions.


        if flag == 0: #MDA
            FileDialog = Gtk.FileChooserDialog(title = "Choose MDA Folder", action = Gtk.FileChooserAction.SELECT_FOLDER)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                self.MDA_folder = FileDialog.get_filename()
                if "mda" in self.MDA_folder:
                    self.Folder_Open()
                    self.Image_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir)),"Images")
                    self.h5_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir)),"h5")                
                    try:
                        self.sheet = self.client.open("S26_"+self.MDA_folder.split(os.sep)[-2])
                    except Exception as e:
                        print(e)
                        self.sheet = None
                    else:
                        self.sheet = self.sheet.get_worksheet(0)
                FileDialog.destroy()
            else:
                FileDialog.destroy()
        elif flag == 1: # previously for selecting Images folder is not found, not useless
            FileDialog = Gtk.FileChooserDialog(title = "Choose Image Folder", action = Gtk.FileChooserAction.SELECT_FOLDER)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                self.Image_folder = FileDialog.get_filename()
            FileDialog.destroy()
        elif flag == 2:
            FileDialog = Gtk.FileChooserDialog(title = "Choose a file", action = Gtk.FileChooserAction.OPEN)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                filename = FileDialog.get_filename()
            else:
                filename = None
            FileDialog.destroy()
            return filename
        elif flag == 3:
            FileDialog = Gtk.FileChooserDialog(title = "Choose a file", action = Gtk.FileChooserAction.SAVE)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                filename = FileDialog.get_filename()
            else:
                filename = None
            FileDialog.destroy()
            return filename
        elif flag == 4:
            FileDialog = Gtk.FileChooserDialog(title = "Choose an HDF5 file", action = Gtk.FileChooserAction.OPEN)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            h5_filter = Gtk.FileFilter()
            h5_filter.set_name("HDF5 files")
            h5_filter.add_pattern("*.h5")
            h5_filter.add_pattern("*.hdf5")
            FileDialog.add_filter(h5_filter)
            all_filter = Gtk.FileFilter()
            all_filter.set_name("All files")
            all_filter.add_pattern("*")
            FileDialog.add_filter(all_filter)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                filename = FileDialog.get_filename()
            else:
                filename = None
            FileDialog.destroy()
            return filename


    def Detector_Changed(self, widget):
        # Updates global detector dimensions and resets display axes when a new detector is selected.


        if not widget.get_active():
            return
        text = widget.get_label()
        self.pilatus_enabled = False
        self.eiger_enabled = False
        if text == "1028 x 1062":
            self.dimY = 1062
            self.dimX = 1028
            self.eiger_enabled = True
        elif text == "516 x 516":
            self.dimY = 516
            self.dimX = 516
        elif text == "515 x 515":
            self.dimY = 515
            self.dimX = 515
        elif text == "256 x 256":
            self.dimY = 256
            self.dimX = 256
        elif text == "487 x 195":
            self.dimY = 195
            self.dimX = 487
            self.pilatus_enabled = True
        elif text == "1360 x 1024":
            self.dimY = 1024
            self.dimX = 1360
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
        self.Image_Axe.cla()
        self.Image_Axe.set_axis_off()
        del self.Image_Image
        #self.Image_Axe.xaxis.set_ticklabels([])
        #self.Image_Axe.yaxis.set_ticklabels([])


    def Change_Colormap(self, widget):
        # Changes the colormap used for all image and 2D plot displays.


        if widget.get_active(): #the callback is invoked two times, deselecting one and selecing another
            self.cm = widget.get_label()
            if not self.xrf_mode:
                self.Image_Image.set_cmap(self.cm)
                self.Image_Canvas.draw()
            self.Plot2D_Image.set_cmap(self.cm)
            self.Plot2D_Canvas.draw()

    def SparseToggled(self, widget):
        # Toggles the use of sparse matrices for efficient memory handling of large image stacks.


        if widget.get_active():
            self.sparse_enabled = True
        else:
            self.sparse_enabled = False

    def ShowAngleToggled(self, widget):
        # Toggles the display of calculated 2-theta/gamma angles instead of pixel coordinates.


        if widget.get_active():
            analysis_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir, 'Analysis')))    
            if os.path.exists(os.path.join(analysis_folder, "gamma.csv")):
                self.gamma = np.genfromtxt(os.path.join(analysis_folder, "gamma.csv"), delimiter=',')
                self.twotheta = np.genfromtxt(os.path.join(analysis_folder, "twotheta.csv"), delimiter=',')
            else:
                self.gamma = np.zeros((1062,1028))
                self.twotheta = np.copy(self.gamma)
            self.show_angle = True
        else:
            self.show_angle = False

    def DirtyFixToggled(self, widget):
        # Toggles a data correction for specific detector types where the first point might be corrupted.

        
        if widget.get_active():
            self.dirty_fix = True
        else:
            self.dirty_fix = False

    def PumpProbeToggled(self, widget):
        # Toggles specialized processing for pump-probe experiments (e.g., interleaving/splitting images).

        
        if widget.get_active():
            self.pump_probe = True
        else:
            self.pump_probe = False

    def XRFModeToggled(self, widget):
        # Switches the UI layout between diffraction image mode and XRF spectrum mode.

        
        if widget.get_active():
            self.xrf_mode = True
            self.Image_Notebook.set_current_page(1)
        else:
            self.xrf_mode = False
            self.Image_Notebook.set_current_page(0)

    def MainWindow_Destroy(self, widget): 
        # Cleanup routine: closes files, destroys sub-windows, and exits the main GTK loop.

        self.image = None
        try:
            self.h5.close()
        except:
            pass
        try:
            self.Shortcuts.win.destroy()
        except:
            pass
        try:
            self.ProbeLibrary.win.destroy()
        except:
            pass
        try:
            self.PowderHelper.win.destroy()
        except:
            pass
        try:
            self.ShiftCorrection.win.destroy()
        except:
            pass
        try:
            self.Analyze5D.win.destroy()
        except:
            pass
        try:
            self.PtychoLib.win.destroy()
        except:
            pass
        try:
            self.ShowMetadata.win.destroy()
        except:
            pass
        self.Main_Window.destroy()
        Gtk.main_quit()


    def AboutThisProgram(self, widget):
        # Displays the "About" dialog with versioning and developer information.


        dialog = Gtk.AboutDialog()
        dialog.set_program_name("dude")
        dialog.set_comments("Diffraction User Data Explorer")
        dialog.set_version("2.7.0")
        dialog.set_website("https://hub.docker.com/r/danielzt12/dude")
        dialog.set_website_label("https://hub.docker.com/r/danielzt12/dude")
        dialog.get_child().get_children()[1].hide()
        dialog.get_child().get_children()[0].get_children()[0].hide()
        dialog.run()
        

    #---------------scan toolbox callbacks--------------------#   
    # Methods for ROI management and hot pixel masking tools

    def Scan_ToolBox_ImageData_MaskAbove(self, widget):
        # Generates a custom pixel mask based on a hot-pixel threshold.

        if widget.get_active():
            threshold = int(self.Scan_ToolBox_ImageData_HotPixel_Adjustment.get_value())
            self.custom_mask = self.Image_Image.get_array()<threshold
        else:
            self.custom_mask = None


    def Scan_ToolBox_ImageData_MaskBelow(self, widget):
        # Zeroes out all pixel values below the specified threshold across all loaded images.

        threshold = int(self.Scan_ToolBox_ImageData_HotPixel_Adjustment.get_value())
        if self.sparse_enabled:
            self.image.data[self.image.data<threshold]=0
            self.image.eliminate_zeros()
        else:
            self.image[self.image<threshold] = 0


    def Scan_ToolBox_U_Toggled(self, widget):
        # Toggles the visualization of the currently selected Region of Interest (ROI).

        if not widget.get_active():
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.Image_Canvas.draw()
        self.Scan_ToolBox_Plot_Changed(None)


    def Scan_ToolBox_CustomROI_Summed(self, widget):
        # Integrates intensity within the defined ROI across all points in the scan for 1D/2D plotting.

        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        xmin = int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        xmax = int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        ymin = int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        ymax = int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        if ndim == 1:
            xdata =  np.array(self.data[1].p[0].data)
            if self.sparse_enabled:
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata = np.array(self.image[:,indice0].sum(1))[:,0]
            else:
                ydata = self.image[:, ymin:ymax+1,xmin:xmax+1].sum(1).sum(1)
            self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata[(xdata!=0)*(ydata!=0)])
            if self.pump_probe:
                if ymin>550:
                    ymin -= 550
                    ymax -= 550
                else:
                    ymin += 550
                    ymax += 550
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata2 = np.array(self.image[:,indice0].sum(1))[:,0]
                self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata2[(xdata!=0)*(ydata!=0)], hold=True) 
            self.Plot1D_Canvas.draw()
        elif ndim == 2:
            xdata2 = np.zeros((self.data[0]["dimensions"][0],self.data[0]["dimensions"][1]))
            xdata1 = np.copy(xdata2)+1
            ydata = np.copy(xdata2)
            datatmp = np.array(self.data[2].p[0].data)
            xdata2[:datatmp.shape[0]] = datatmp
            xdata1 *= np.array(self.data[1].p[0].data)[:,np.newaxis]
            if self.sparse_enabled:
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata = self.image[:,indice0].sum(1).reshape(xdata1.shape)
            else:
                ydata = self.image[:, ymin:ymax+1,xmin:xmax+1].sum(1).sum(1).reshape(xdata1.shape)
            self.Plot2D_Axe.cla()
            self.Plot2D_Axe.set_axis_off()
            if self.pump_probe:
                if ymin>550:
                    ymin -= 550
                    ymax -= 550
                else:
                    ymin += 550
                    ymax += 550
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                # for some reasons the images aren't right, and needs to be shifted
                ydata2 = self.image[:,indice0].sum(1).reshape(xdata1.shape)
                ydata_flat = ydata.flatten()
                ydata_flat[0,:-1] = ydata_flat[0,1:]
                ydata = ydata_flat.reshape(ydata.shape)
                ydata_flat = ydata2.flatten()
                ydata_flat[0,:-1] = ydata_flat[0,1:]
                ydata2 = ydata_flat.reshape(ydata2.shape)
                ydata = np.hstack((ydata,ydata2))
                self.Plot2D_xdata1 = np.hstack((self.Plot2D_xdata1,self.Plot2D_xdata1))
                self.Plot2D_xdata2 = np.hstack((self.Plot2D_xdata2,self.Plot2D_xdata2))
            self.Plot2D_Image = self.Plot2D_Axe.imshow(ydata,  aspect = 'equal', interpolation = 'nearest', cmap=self.cm, origin="lower")
            self.Plot2D_ydata = ydata
            if self.Plot2D_Log_ToggleButton.get_active():
                self.Plot2D_Image.set_norm(colors.LogNorm())
            else:
                self.Plot2D_Image.set_norm(colors.Normalize())
            self.Plot2D_Canvas.draw()
        
    def Scan_ToolBox_CustomROI_Added(self, widget):
        # Saves the current ROI coordinates into the custom ROI list for later reuse.

        
        roi_name = self.Scan_ToolBox_CustomROI_Entry.get_text()
        roi_info = [roi_name, \
                        int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())]
        for row in self.CustomROI_store:
            if roi_info[0] == row[0]:
                roi_info[0]+="*"
        self.CustomROI_store.append(roi_info)
        self.Scan_ToolBox_U_ComboBox.set_active(len(self.Scan_ToolBox_U_ComboBox.get_model())-1)
        self.Scan_ToolBox_U_ToggleButton.set_active(True)


    def Scan_ToolBox_CustomROI_Changed(self, widget):
        # Updates the ROI coordinate limits and visual overlay when spin buttons are adjusted.

        
        xmin = int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        xmax = int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        ymin = int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        ymax = int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(xmax-1)
        self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(xmin+1)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(ymax-1)
        self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(ymin+1)
        
        try:
            self.Image_ShowROI_Rectangle.remove()
        except (AttributeError, ValueError):
            pass
        self.Image_ShowROI_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b')
        self.Image_Canvas.draw()
        
     
    def Scan_ToolBox_Plot_Changed(self, widget):
        # Refreshes the 1D/2D plots based on the current detector, monitor, or ROI selection.


        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        if ndim == 1:
            xdata =  np.array(self.data[1].p[0].data)

            if self.Scan_ToolBox_U_ToggleButton.get_active():
                roi_info = self.CustomROI_store[self.Scan_ToolBox_U_ComboBox.get_active()]
                if self.sparse_enabled:   # disabling sparse for 1d
                    indice0 = np.zeros((1062,1028))
                    indice0[roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1]=1
                    indice0 = indice0.flatten().nonzero()[0]
                    ydata = np.array(self.image[:,indice0].sum(1))[:,0]
                    #print(ydata.shape)
                else:
                    ydata = self.image[:, roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1].sum(1).sum(1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(roi_info[1])
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(roi_info[2])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(roi_info[3])
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(roi_info[4])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
            else:
                ydata =  np.array(self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].data)
                # dirty fix first point bug
                if self.dirty_fix and \
                   ("eiger" in self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name or \
                   "QMPX3" in self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name ):
                    ydata[:-1] = ydata[1:]
                    # print("dirty fixing first point bug")
            if self.Scan_ToolBox_M_ToggleButton.get_active():
                mdata = np.array(self.data[1].d[self.Scan_ToolBox_M_ComboBox.get_active()].data)
                mdata /= mdata.mean()
                ydata /= mdata
            self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata[(xdata!=0)*(ydata!=0)]) 
            self.Plot1D_Canvas.draw()
        else:
            # this is to avoid crashing on unfinished data
            xdata2 = np.zeros((self.data[0]["dimensions"][0],self.data[0]["dimensions"][1]))
            xdata1 = np.copy(xdata2)+1
            ydata = np.copy(xdata2)
            mdata = np.copy(xdata2)
            datatmp = np.array(self.data[2].p[0].data)
            xdata2[:datatmp.shape[0]] = datatmp
            xdata1 *= np.array(self.data[1].p[0].data)[:,np.newaxis]

            if self.Scan_ToolBox_U_ToggleButton.get_active():
                roi_info = self.CustomROI_store[self.Scan_ToolBox_U_ComboBox.get_active()]
                if self.sparse_enabled:
                    indice0 = np.zeros((1062,1028))
                    indice0[roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1]=1
                    indice0 = indice0.flatten().nonzero()[0]
                    ydata = self.image[:,indice0].sum(1).reshape(xdata1.shape)
                else:
                    ydata = self.image[:, roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1].sum(1).sum(1).reshape(xdata1.shape)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(roi_info[1])
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(roi_info[2])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(roi_info[3])
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(roi_info[4])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
            else:
                datatmp =  np.array(self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].data)
                ydata[:datatmp.shape[0]] = datatmp
                if self.dirty_fix and \
                   ("eiger" in self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name or \
                    "QMPX3" in self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name ) and \
                   self.data[1].time != "whatever":
                    print("dirty fixing first point bug")
                    ydata_flat = ydata.flatten()
                    ydata_flat[:-1] = ydata_flat[1:]
                    ydata = ydata_flat.reshape(ydata.shape)
            if self.Scan_ToolBox_M_ToggleButton.get_active():
                datatmp = np.array(self.data[2].d[self.Scan_ToolBox_M_ComboBox.get_active()].data)
                mdata[:datatmp.shape[0]] = datatmp
                mdata /= mdata[mdata>0].mean()
                mdata[mdata==0] = 1
                ydata /= mdata

            ydata_min = ydata[ydata!=0].min()
            self.Plot2D_Vmin_HScale_Adjustment.set_lower(ydata_min)
            self.Plot2D_Vmax_HScale_Adjustment.set_upper(ydata.max())
            if self.Plot2D_AutoScale_ToggleButton.get_active():
                self.Plot2D_Vmin_HScale_Adjustment.handler_block(self.Plot2D_Vmin_Changed_Handler)
                self.Plot2D_Vmax_HScale_Adjustment.handler_block(self.Plot2D_Vmax_Changed_Handler)
                self.Plot2D_Vmin_HScale_Adjustment.set_upper(ydata.max()-.01)
                self.Plot2D_Vmin_HScale_Adjustment.set_value(ydata_min)
                self.Plot2D_Vmax_HScale_Adjustment.set_lower(ydata_min+.01)
                self.Plot2D_Vmax_HScale_Adjustment.set_value(ydata.max())
                self.Plot2D_Vmin_HScale_Adjustment.handler_unblock(self.Plot2D_Vmin_Changed_Handler)
                self.Plot2D_Vmax_HScale_Adjustment.handler_unblock(self.Plot2D_Vmax_Changed_Handler)
            if self.MDA_File_ListStore[self.mda_selection_path[0]][1] != "ao01":
                self.Scan_Plot2D(xdata1, xdata2, ydata)
            else:
                self.Scan_PlotSpiral(xdata1, xdata2, ydata)
            self.Plot2D_Canvas.draw()


    def Scan_TreeView_Selection_Changed(self, treeselection):
        # Triggered when a different MDA file is selected in the file list.


        treemodelfilter, filterpathlist = treeselection.get_selected_rows()
        self.mda_selection_path = list(map(treemodelfilter.convert_path_to_child_path, filterpathlist))
        mdapath = self.MDA_File_ListStore[self.mda_selection_path[0]][-2]
        self.Scan_ToolBox_U_ToggleButton.set_active(False)
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(False)
        self.AllImagesLoaded = False
        self.Scan_Load(mdapath)


    def Scan_TreeView_Select_Hijack(self, treeselection, treefilter, irow, selected):
        # Custom selection logic to allow browsing files without automatically triggering heavy data loads.

        
        irow, icolumn = self.MDA_File_TreeView.get_cursor()
        if icolumn == self.MDA_File_TreeView.get_column(13): 
            self.MDA_File_ListStore[irow][-1] = not self.MDA_File_ListStore[irow][-1]
            return False
        else:
            return True


    def Scan_Load(self, mdapath):
        # Core logic for loading a scan's metadata, 1D/2D data, and associated detector image files.
        # Handles complex cases like interrupted scans, missing images, and various detector formats.

        Scan_ToolBox_Y_ComboBox_Select = self.Scan_ToolBox_Y_ComboBox.get_active()
        Scan_ToolBox_M_ComboBox_Select = self.Scan_ToolBox_M_ComboBox.get_active()
        self.Scan_ToolBox_Y_ComboBox.handler_block(self.Scan_ToolBox_Y_Changed_Handler)
        self.Scan_ToolBox_M_ComboBox.handler_block(self.Scan_ToolBox_M_Changed_Handler)
        self.MDA_Det_store.clear()
        self.data = readMDA(mdapath, verbose=0, maxdim=3 if self.xrf_mode else 2)
        self.xrf_data = []
        xrf_in_mda = self.data[0]["dimensions"][-1] == 2048
        ndim = len(self.data)-2 if xrf_in_mda else len(self.data)-1
        xrf_in_netcdf = False
        if self.xrf_mode:
            if xrf_in_mda:
                nchan = self.data[ndim+1].nd
                for i in range(nchan):
                    self.xrf_data += [self.data[ndim+1].d[i].data]
                self.xrf_data = np.array(self.xrf_data).reshape(nchan,-1,2048)
            else:
                fluo_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir, 'fluo')))    
                f_fluo = [f for f in os.listdir(fluo_folder) if f.startswith("scan_{0}".format(self.MDA_File_ListStore[self.mda_selection_path[0]][0]))]
                if len(f_fluo):
                    f_fluo = f_fluo[0]
                    netcdffile = netCDF4.Dataset(os.path.join(fluo_folder,f_fluo), "r")
                    ny = self.data[0]["dimensions"][0]
                    nx = self.data[0]["dimensions"][1]
                    self.xrf_data = np.swapaxes(np.swapaxes((netcdffile.variables['array_data'][:,:,256:].reshape(ny, 2, 124, 256+2048*4)[:,:,:nx,256:]),1,2).reshape(-1,8,2048), 0,1)
                    netcdffile.close()
                    xrf_in_netcdf = True
        for i in range(self.data[ndim].nd):
            dname = self.data[ndim].d[i].name
            self.MDA_Det_store.append([dname.replace("s26_eiger_cnm","eiger"), 1])
        self.Scan_ToolBox_M_ComboBox.set_active(Scan_ToolBox_M_ComboBox_Select)
        self.Scan_ToolBox_Y_ComboBox.handler_unblock(self.Scan_ToolBox_Y_Changed_Handler)
        self.Scan_ToolBox_M_ComboBox.handler_unblock(self.Scan_ToolBox_M_Changed_Handler)
        self.Plot_Notebook.set_current_page(ndim-1)        
        self.Scan_ToolBox_Y_ComboBox.set_active(Scan_ToolBox_Y_ComboBox_Select)

        if self.xrf_mode and (xrf_in_mda or xrf_in_netcdf):
            nspec = self.xrf_data.shape[1] if ndim == 1 else self.xrf_data.shape[1]*self.xrf_data.shape[2]
            self.XRF_Plot_HScale_Adjustment.set_upper(nspec-1)
            self.XRF_Plot_HScale_Adjustment.set_lower(0)
            self.XRF_Plot_HScale_Adjustment.set_value(0)
        else:
            if self.eiger_enabled:
                h5_filename = [f.name for f in os.scandir(self.h5_folder) if "scan_{0}".format(self.MDA_File_ListStore[self.mda_selection_path[0]][0]) in f.name]
                if len(h5_filename) == 1:
                    h5_filename = h5_filename[0]
                else:
                    h5_index = np.array([int(f.split(".")[0].split("_")[-1]) for f in h5_filename])
                    h5_filename = np.take(h5_filename, h5_index.argsort())[0]
                try:
                    self.h5.close()
                except Exception as e:
                    pass
                try:
                    self.h5 = h5py.File(os.path.join(self.h5_folder, h5_filename), 'r')
                except Exception as e:
                    self.image_index_max = -1
                else:
                    self.image_index_max = self.h5['entry/data/data'].shape[0]-1
                    self.image_index_min = 0
                    self.nbin = 1
                    if hasattr(self, "ShowMetadata"):
                        self.ShowMetadata.set_h5(self.h5)
                    # Extract coarse motor positions for display
                    self._update_coarse_motors()

            else:
                self.image_path = os.path.join(self.Image_folder, mdapath.split(".")[0].split("_")[-1].lstrip("0"))
                # os.chmod(self.image_path, 0o777)
                if self.pilatus_enabled == True:
                    self.image_list = [imagefile.name for imagefile in os.scandir(self.image_path) if ('.tif' in imagefile.name and 'il' in imagefile.name)]
                else:
                    self.image_list = [imagefile.name for imagefile in os.scandir(self.image_path) if ('.tif' in imagefile.name and 'il' not in imagefile.name)]
                self.tif_index = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in self.image_list])
                self.image_list = np.take(self.image_list , self.tif_index.argsort())
                self.nbin = 1
                if ndim == 2:
                    self.ntiff = self.data[0]['dimensions'][0] * self.data[0]['dimensions'][1]
                    for i in range(self.data[ndim].nd):
                        dname = self.data[ndim].d[i].name
                        if "FileNumber" in dname:
                            mda_index = np.array(self.data[2].d[i].data).flatten()
                    #print("a,",mda_index)
                    print (mda_index.max(), self.tif_index.max(), mda_index.min(), self.tif_index.min())
                else:
                    self.ntiff = self.data[0]['dimensions'][0]
                if len(self.image_list) != self.ntiff:
                    print ("I only print this because the file saving here has been messy!")
                    print ("Expecting {0} files but got {1} instead".format(self.ntiff, len(self.image_list)))
                    if ndim == 1:
                        if self.ntiff > len(self.image_list):
                            print ("I suppose the scan is interrupted and will not attempt to correct that")
                        else:
                            if len(self.image_list)%self.ntiff == 0:
                                self.nbin = max(1,int(len(self.image_list)/self.ntiff))
                                print ("I suppose the scan has {0} repeated measurement for each point".format(self.nbin))
                                self.image_list = self.image_list[0::self.nbin]
                            else:
                                print ("There is a chance that files are being saved into the wrong folder. Please correct manually 1!")
                    else:
                        if len(self.image_list) > self.ntiff:
                            self.nbin = max(1,int(len(self.image_list)/self.ntiff))
                            if len(self.image_list)%self.ntiff:
                                self.nbin+=1
                            if self.nbin>1:
                                if not len(self.image_list)%self.ntiff:
                                    print ("I suppose the scan has {0} repeated measurement for each point".format(self.nbin))
                                    self.image_list = self.image_list[0::self.nbin]
                                else:
                                    if mda_index.min() == self.tif_index.min():
                                        inc_index = mda_index[1:] - mda_index[:-1]
                                        hiccup_index = (inc_index == 2).nonzero()[0]
                                        if inc_index.max() == 2:
                                            print ("I suppose the scan has {0} repeated measurement for each point, but there are {1} images missing here!".format(self.nbin, self.ntiff*self.nbin-len(self.image_list)))
                                            if inc_index[hiccup_index+1].sum() == 0 : # because a hiccup is always followed by a 0                                 
                                                print("Experiencing {0} hiccups. Use only tif index.".format((inc_index==2).sum()))
                                                mda_index[hiccup_index+1] -= 1
                                                self.image_list = self.image_list[(mda_index-mda_index.min()).astype(int)]
                                                self.image_list = self.image_list[0::self.nbin]
                                            else:
                                                print("Experiencing {0} hiccups, unable to correct.".format((inc_index==2).sum()))
                                        elif inc_index.min() < 0 or inc_index.max() > 2:
                                            print ("This should not happen. Please correct manually 5!")
                                        elif inc_index.min() == 0 and (inc_index==0).sum() == mda_index.shape[0] - self.tif_index.shape[0] + hiccup_index.shape[0]:
                                            print("Lost {0} images\nDeciding to trust mda index.".format((inc_index==0).sum()- hiccup_index.shape[0]))
                                            self.image_list = self.image_list[0::self.nbin]
                                        else:
                                            print("Dunno what to do 7!")
                                    else:
                                        print("Dunno what to do 6!")
                            else:
                                if mda_index.max() < self.tif_index.max():
                                    if mda_index.min() == self.tif_index.min():
                                        print ("Ignoring {0} extra images in the folder.\nThey might belong to the next scan.".format(self.tif_index.max() - mda_index.max()))
                                        self.image_list = self.image_list[self.tif_index<=mda_index.max()]
                                        self.tif_index = self.tif_index[self.tif_index<=mda_index.max()]
                                    else:
                                        print ("There is a chance that files are being saved into the wrong folder. Please correct manually 2!")
                                else:
                                    if mda_index.max() == self.tif_index.max() and mda_index.min() == self.tif_index.min() +1:
                                        print ("This bug is new in 2019, the mda index is shifted by 1 with regard to the tif index, but strangely for these long scans there is one more point in the beginning")
                                        self.image_list = self.image_list[:-1]
                                    else:
                                        print ("This should not happen. Please correct manually 3!")
                        else:
                            if mda_index.min() == self.tif_index.min():
                                inc_index = mda_index[1:] - mda_index[:-1]
                                hiccup_index = (inc_index == 2).nonzero()[0]
                                if inc_index.min() == 0 and (inc_index==0).sum() == mda_index.shape[0] - self.tif_index.shape[0] + hiccup_index.shape[0]:
                                    print("Lost {0} images\nDeciding to trust mda index.".format((inc_index==0).sum()- hiccup_index.shape[0]))
                                    self.image_list = self.image_list[(mda_index-mda_index.min()).astype(int)]
                                elif inc_index.max() == 2:
                                    if inc_index[hiccup_index+1].sum() == 0 : # because a hiccup is always followed by a 0
                                        print("Experiencing {0} hiccups. Use only tif index.".format((inc_index==2).sum()))
                                        mda_index[hiccup_index+1] -= 1
                                        self.image_list = self.image_list[mda_index]
                                    else:
                                        print("Experiencing {0} hiccups, unable to correct.".format((inc_index==2).sum()))
                                elif inc_index.min() < 0 or inc_index.max() > 2:
                                    print ("This should not happen. Please correct manually 4!")
                                else:
                                    print("Dunno what to do 8!")
                            else:
                                if mda_index.min() == self.tif_index.min() +1:
                                    print ("This bug is new in 2019, the mda index is shifted by 1 with regard to the tif index, and I suppose you interrupted the scan. Other than that all is well")
                                else:
                                    print ("There is a chance that files are being saved into the wrong folder. Please correct manually!")
                if self.tif_index.shape[0]:
                    self.image_index_min = self.tif_index.min()
                    self.image_index_max = self.tif_index.max()
                    ###self.image_index_min = mda_index.min()
                    ###self.image_index_max = mda_index.max()

            self.Image_Plot_HScale_Adjustment.set_upper(self.image_index_max)
            self.Image_Plot_HScale_Adjustment.set_lower(self.image_index_min)
            if self.Image_Plot_HScale_Adjustment.get_value() == 0:
                self.Image_Plot_HScale_Changed(self.Image_Plot_HScale_Adjustment)
            else:
                self.Image_Plot_HScale_Adjustment.set_value(self.image_index_min)
            self.Image_Plot_HScale_Adjustment.set_step_increment(self.nbin)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(False)
        

    def _update_coarse_motors(self):
        """Read SAMX/SAMY from the open H5 file and update the display label + table row."""
        sx_str, sy_str = '---', '---'
        try:
            sx = self.h5['/entry/instrument/26-ID-C/SAMX'][()]
            sx_str = '{0:.1f}'.format(float(sx.flat[0]) if hasattr(sx, 'flat') else float(sx))
        except Exception:
            pass
        try:
            sy = self.h5['/entry/instrument/26-ID-C/SAMY'][()]
            sy_str = '{0:.1f}'.format(float(sy.flat[0]) if hasattr(sy, 'flat') else float(sy))
        except Exception:
            pass
        sz_str = '---'
        try:
            sz = self.h5['/entry/instrument/26-ID-C/SAMZ'][()]
            sz_str = '{0:.1f}'.format(float(sz.flat[0]) if hasattr(sz, 'flat') else float(sz))
        except Exception:
            pass
        # Update the label below the 2D plot
        self.Plot2D_CoarseMotor_Label.set_markup(
            '<span font_desc="monospace 9" foreground="#555555">'
            'SAMX: {0}   SAMY: {1}   SAMZ: {2}</span>'.format(sx_str, sy_str, sz_str))
        # Also update the scan table columns for this row
        try:
            row = self.mda_selection_path[0]
            self.MDA_File_ListStore[row][10] = sx_str
            self.MDA_File_ListStore[row][11] = sy_str
            self.MDA_File_ListStore[row][12] = sz_str
        except Exception:
            pass

    def Folder_Refresh(self, widget):
        # Re-scans the MDA folder for new files added since the last update.

        self.MDA_cursor_current = self.MDA_File_TreeView.get_cursor()[0] #unstable
        self.MDA_File_TreeView.get_selection().handler_block(self.MDA_File_TreeView_Selection_Changed_Handler) #unstable
        MDAfile_list = [mdafile for mdafile in os.listdir(self.MDA_folder) if mdafile.endswith('.mda')]
        index_list = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in MDAfile_list])
        last_index = int(self.MDA_File_ListStore[-1][0])
        MDAfile_list = np.array(MDAfile_list)[index_list>=last_index]
        index_list = index_list[index_list>=last_index]
        MDAfile_list = np.take(MDAfile_list , index_list.argsort())
        self.Folder_Scan(MDAfile_list) 


    def Folder_Open(self):
        # Clears the current file list and performs a full scan of the MDA directory.

        self.MDA_File_TreeView.get_selection().handler_block(self.MDA_File_TreeView_Selection_Changed_Handler)
        self.MDA_File_ListStore.clear()
        self.MDA_File_TreeView.get_selection().handler_unblock(self.MDA_File_TreeView_Selection_Changed_Handler)
        MDAfile_list = [mdafile for mdafile in os.listdir(self.MDA_folder) if mdafile.endswith('mda')]
        index_list = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in MDAfile_list])
        MDAfile_list = np.take(MDAfile_list , index_list.argsort())
        self.Folder_Scan(MDAfile_list)


    def Folder_Scan(self, MDAfile_list):
        # Parallelized loading of MDA metadata for a list of files to populate the tree view.

        if len(self.MDA_File_ListStore):
            self.MDA_File_ListStore.remove(self.MDA_File_ListStore[-1].iter)
        mda_list = list(map(lambda x:os.path.join(self.MDA_folder,x), MDAfile_list))
        pool = Pool(processes=cpu_count()) 
        for i, tmplist in enumerate(pool.imap(mda_loader, mda_list)):
            self.MDA_File_ListStore.append(tmplist)
        pool.close()
        self.MDA_File_TreeView.get_selection().handler_unblock(self.MDA_File_TreeView_Selection_Changed_Handler)#unstable
        if hasattr(self, "MDA_cursor_current"):
            self.MDA_File_TreeView.set_cursor(self.MDA_cursor_current)#unstable


    def Upload_To_Logbook(self, widget, flag):
        # Saves the current active plot as a JPEG and uploads it to Google Drive/Logbook.

        folder_id = '1HaYq1NBpA5CTgdP1Y__rOPwpbRoeL6Vo'
        filename = str(uuid.uuid4())+'.jpg'
        file_metadata = {
            'name': filename,
            'parents': [folder_id]
        }
        filename = os.path.join("/tmp", filename)
        if flag == 0:
            self.Image_Figure.savefig(filename, bbox_inches='tight')
            self.Image_Canvas.draw() # somehow the image will just disappear, try to fix that
        elif flag == 1:
            self.Plot1D_Figure.savefig(filename, bbox_inches='tight')
            self.Plot1D_Canvas.draw()  # somehow the image will just disappear, try to fix that
        elif flag == 2:
            self.Plot2D_Figure.savefig(filename, bbox_inches='tight')
            self.Plot2D_Canvas.draw()  # somehow the image will just disappear, try to fix that
        else:
            self.XRF_Figure.savefig(filename, bbox_inches='tight')
            self.XRF_Canvas.draw()  # somehow the image will just disappear, try to fix that
        media = MediaFileUpload(filename,
                                mimetype='image/jpeg',
                                resumable=False)
        newfile = self.service.files().create(body=file_metadata,
                                              media_body=media,
                                              fields='id').execute()
                    
        fileid = newfile.get('id')
        permission = {'role': 'reader', 'type': 'anyone'}
        self.service.permissions().create(fileId=fileid,body=permission).execute()
        weblink = self.service.files().get(fileId=fileid,fields='webContentLink').execute()
        scannum = self.MDA_File_ListStore[self.mda_selection_path[0]][0]
        i_row = self.sheet.col_values(3).index(scannum)+2
        #print(i_row)
        content = "=IMAGE(\"{0}\", 1)".format(weblink['webContentLink'].replace("'",'"'))
        self.sheet.insert_row(["","",str(scannum),"",content], i_row, value_input_option='USER_ENTERED')
        set_row_height(self.sheet, str(i_row), 200)
               

    def Console_KeyPressed(self, widget, event):
        # Implements history navigation (Up/Down arrows) for the GUI python console.

        keyname = Gdk.keyval_name(event.keyval)
        if keyname == 'Up' and (event.state & Gdk.ModifierType.SHIFT_MASK):
            if self.console_cmd_index > 0:
                self.console_cmd_index -= 1
            else:
                self.console_cmd_index = len(widget.get_completion().get_model())-1
            widget.set_text(widget.get_completion().get_model()[self.console_cmd_index][0])
        elif keyname == 'Down' and (event.state & Gdk.ModifierType.SHIFT_MASK):
            if self.console_cmd_index < len(widget.get_completion().get_model()) -1:
                self.console_cmd_index += 1
            else:
                self.console_cmd_index = 0
            widget.set_text(widget.get_completion().get_model()[self.console_cmd_index][0])
 

    def __init__(self):
        # Constructor for the main window. Initializes state, UI components, and authentication.
        # Authenticate with Google Services for logbook integration if configuration exists
        if os.path.exists("/home/sector26/pythonscripts/Tao/token.json"):
            SCOPES = ['https://www.googleapis.com/auth/drive.file']
            creds = Credentials.from_authorized_user_file('/home/sector26/pythonscripts/Tao/token.json', SCOPES)
            # If there are no (valid) credentials available, let the user log in.
            if not creds or not creds.valid:
                if creds and creds.expired and creds.refresh_token:
                    creds.refresh(Request())
                else:
                    flow = InstalledAppFlow.from_client_secrets_file('/home/sector26/pythonscripts/Tao/client_secrets.json', SCOPES)
                    creds = flow.run_local_server()
                    # Save the credentials for the next run
                with open('/home/sector26/pythonscripts/Tao/token.json', 'w') as token:
                    token.write(creds.to_json())

            self.service = build('drive', 'v3', credentials=creds)

            SCOPES2 = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']

            creds2 = ServiceAccountCredentials.from_json_keyfile_name('/home/sector26/pythonscripts/Tao/automatedlogbook.json', SCOPES2)

            self.client = gspread.authorize(creds2)
            beamline = True
        else:
            beamline = False

        # Initialize default state variables
        self.cm = "viridis"

        self.dimY = 1062
        self.dimX = 1028
        self.Image_Rectangle_Drawing = False
        self.Image_Zoomed = False
        self.Image_Zooming = False
        self.Image_Panning = False
        self.non_animated_background = None
        self.custom_mask = None

        #------------------------MDA File Window------------------------------#    
        # Setup the scrolled window and list store for the MDA file browser

        MDA_File_ScrolledWindow = Gtk.ScrolledWindow()
        MDA_File_ScrolledWindow.set_policy(Gtk.PolicyType.ALWAYS, Gtk.PolicyType.ALWAYS)
        MDA_File_ScrolledWindow.set_overlay_scrolling(False)
        self.MDA_File_ListStore = Gtk.ListStore(str, str, str, str, str, str, str, str, str, str, str, str, str, str, bool)
        self.MDA_File_TreeView_Filter = self.MDA_File_ListStore.filter_new()
        self.MDA_File_TreeView = Gtk.TreeView.new_with_model(self.MDA_File_TreeView_Filter)
        self.MDA_File_TreeView.get_selection().set_mode(Gtk.SelectionMode.SINGLE)
        self.MDA_File_TreeView_Selection_Changed_Handler = self.MDA_File_TreeView.get_selection().connect("changed", self.Scan_TreeView_Selection_Changed)
        # hijack the default function to allow multiselect without loading all the scans
        self.MDA_File_TreeView.get_selection().set_select_function(self.Scan_TreeView_Select_Hijack)

        i = 0      
        for name in ["N#","M1","MIN","MAX","NP","M2","MIN","MAX","NP","CT","SX","SY","SZ"]:
            MDA_File_CellRendererText = Gtk.CellRendererText()     
            MDA_File_CellRendererText.set_property("xalign", 1)
            MDA_File_TreeViewColumn = Gtk.TreeViewColumn(name, MDA_File_CellRendererText, text = i)
            i += 1
            self.MDA_File_TreeView.append_column(MDA_File_TreeViewColumn)
        MDA_File_CellRendererToggle = Gtk.CellRendererToggle()     
        MDA_File_CellRendererToggle.set_property("xalign", .5)
        MDA_File_TreeViewColumn = Gtk.TreeViewColumn("", MDA_File_CellRendererToggle, active = 14)
        self.MDA_File_TreeView.append_column(MDA_File_TreeViewColumn)
        self.MDA_File_TreeView.set_enable_search(False)
        MDA_File_ScrolledWindow.add(self.MDA_File_TreeView)
        MDA_File_ScrolledWindow.set_size_request(360, 320)

        
        # ROI Plot : X, Y, Mon
        # Setup widgets for selecting Detectors, ROIs, and Monitors for 1D/2D plotting

        self.MDA_Det_store = Gtk.ListStore(str, int)
        renderer_text = Gtk.CellRendererText()    
        self.CustomROI_store = Gtk.ListStore(str, int, int, int, int)
        renderer_text2 = Gtk.CellRendererText()    

        self.Scan_ToolBox_M_ToggleButton = Gtk.ToggleButton(label = " MON ")
        self.Scan_ToolBox_M_ComboBox = Gtk.ComboBox.new_with_model(self.MDA_Det_store)
        self.Scan_ToolBox_M_ComboBox.pack_start(renderer_text, True)
        self.Scan_ToolBox_M_ComboBox.add_attribute(renderer_text, "text", 0)
        self.Scan_ToolBox_M_ComboBox.add_attribute(renderer_text, "visible", 1)
        self.Scan_ToolBox_M_ToggleButton.connect("toggled", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_Y_ToggleButton = Gtk.ToggleButton(label = " DET ")
        self.Scan_ToolBox_Y_ComboBox = Gtk.ComboBox.new_with_model(self.MDA_Det_store)
        self.Scan_ToolBox_Y_ComboBox.pack_start(renderer_text, True)
        self.Scan_ToolBox_Y_ComboBox.add_attribute(renderer_text, "text", 0)
        self.Scan_ToolBox_Y_ComboBox.add_attribute(renderer_text, "visible", 1)
        self.Scan_ToolBox_U_ToggleButton = Gtk.ToggleButton(label = " ROI ")
        self.Scan_ToolBox_U_ComboBox = Gtk.ComboBox.new_with_model(self.CustomROI_store)
        self.Scan_ToolBox_U_ComboBox.pack_start(renderer_text2, True)
        self.Scan_ToolBox_U_ComboBox.add_attribute(renderer_text2, "text", 0)
        self.Scan_ToolBox_U_ToggleButton.connect("toggled", self.Scan_ToolBox_U_Toggled)
        self.Scan_ToolBox_U_ToggleButton.set_active(False)
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(False)

        Scan_Det_HBox1 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox1.set_border_width(3)
        Scan_Det_HBox1.pack_start(self.Scan_ToolBox_Y_ToggleButton, False, False, 3)
        Scan_Det_HBox1.pack_start(self.Scan_ToolBox_Y_ComboBox, True, True, 3)
        Scan_Det_HBox2 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox2.set_border_width(3)
        Scan_Det_HBox2.pack_start(self.Scan_ToolBox_U_ToggleButton, False, False, 3)
        Scan_Det_HBox2.pack_start(self.Scan_ToolBox_U_ComboBox, True, True, 3)
        Scan_Det_HBox3 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox3.set_border_width(3)
        Scan_Det_HBox3.pack_start(self.Scan_ToolBox_M_ToggleButton, False, False, 3)
        Scan_Det_HBox3.pack_start(self.Scan_ToolBox_M_ComboBox, True, True, 3)
        
        self.Scan_ToolBox_Y_Changed_Handler = self.Scan_ToolBox_Y_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_M_Changed_Handler = self.Scan_ToolBox_M_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_U_Changed_Handler = self.Scan_ToolBox_U_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        
         # ROI Management
         # Setup widgets for naming, adding, and summing custom ROIs

        self.Scan_ToolBox_CustomROI_Entry = Gtk.Entry()
        self.Scan_ToolBox_CustomROI_Entry.set_size_request(40, 20)
        self.Scan_ToolBox_CustomROI_Sum_Button = Gtk.Button(label = " Sum ")
        self.Scan_ToolBox_CustomROI_Sum_Button.connect("clicked", self.Scan_ToolBox_CustomROI_Summed)
        self.Scan_ToolBox_CustomROI_Add_Button = Gtk.Button(label = " Add ")
        self.Scan_ToolBox_CustomROI_Add_Button.connect("clicked", self.Scan_ToolBox_CustomROI_Added)

        Scan_ToolBox_ImageData_LoadAll_Button = Gtk.Button(label = 'Load All')
        Scan_ToolBox_ImageData_LoadAll_Button.connect("clicked", self.LoadAllImages)

        Scan_ToolBox_CustomROI_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Scan_ToolBox_CustomROI_HBox.set_border_width(3)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Entry, True, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Add_Button, False, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Sum_Button, False, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(Scan_ToolBox_ImageData_LoadAll_Button, False, True, 0)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(False)
        
        # ROI Definition
        # Setup spin buttons for precise coordinate definition of ROIs

        self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment = Gtk.Adjustment(value = 200, lower = 0, upper = self.dimY-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment = Gtk.Adjustment(value = 100, lower = 0, upper = self.dimY-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment = Gtk.Adjustment(value = 200, lower = 0, upper = self.dimX-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment = Gtk.Adjustment(value = 100, lower = 0, upper = self.dimX-1, step_increment = 1, page_increment = 2, page_size = 0)
        Scan_ToolBox_CustomROI_YMax_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_YMax_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment)
        Scan_ToolBox_CustomROI_YMax_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_YMin_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_YMin_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment)
        Scan_ToolBox_CustomROI_YMin_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_XMax_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_XMax_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment)
        Scan_ToolBox_CustomROI_XMax_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_XMin_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_XMin_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment)
        Scan_ToolBox_CustomROI_XMin_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_Cut_Button = Gtk.Button(label = ' Cut ')

        self.Scan_ToolBox_CustomROI_YMax_SpinButton_Changed = self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed = self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed = self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed = self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)

        Scan_ToolBox_CustomROI_Table = Gtk.Table(n_rows = 3, n_columns = 3, homogeneous = False)
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_YMax_SpinButton, left_attach = 1, right_attach = 2, top_attach = 2, bottom_attach = 3)
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_YMin_SpinButton, left_attach = 1, right_attach = 2, top_attach = 0, bottom_attach = 1)  
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_XMin_SpinButton, left_attach = 0, right_attach = 1, top_attach = 1, bottom_attach = 2)  
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_XMax_SpinButton, left_attach = 2, right_attach = 3, top_attach = 1, bottom_attach = 2)

        # Customized HotPixel Mask
        # Setup threshold controls and masking buttons for detector data cleanup

        self.Scan_ToolBox_ImageData_HotPixel_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 65536, step_increment = 1, page_increment = 10, page_size = 0)
        Scan_ToolBox_ImageData_HotPixel_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_ImageData_HotPixel_SpinButton.set_adjustment(self.Scan_ToolBox_ImageData_HotPixel_Adjustment)
        Scan_ToolBox_ImageData_HotPixel_SpinButton.set_digits(0)
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton = Gtk.ToggleButton(label = 'Mask Above')
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.connect('clicked', self.Scan_ToolBox_ImageData_MaskAbove)
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.set_tooltip_text('When this button is pressed down, all pixels on the current image which are hotter than the defined threshold will be masked, the mask will then be applied to all images as long as this button is kept pressed down.\n This is independent of the static hot pixel correction in Menu->Settings') 

        self.Scan_ToolBox_ImageData_MaskBelow_Button = Gtk.Button(label = 'Mask Below')
        self.Scan_ToolBox_ImageData_MaskBelow_Button.connect('clicked', self.Scan_ToolBox_ImageData_MaskBelow)
        self.Scan_ToolBox_ImageData_MaskBelow_Button.set_tooltip_text('When this button is pressed, all loaded data with value below the set value is set to 0') 

        Scan_ToolBox_ImageData_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Scan_ToolBox_ImageData_HBox.set_border_width(3)
        Scan_ToolBox_ImageData_HBox.pack_start(Scan_ToolBox_ImageData_HotPixel_SpinButton, True, True, 0)
        Scan_ToolBox_ImageData_HBox.pack_start(self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton, True, False, 0)
        Scan_ToolBox_ImageData_HBox.pack_start(self.Scan_ToolBox_ImageData_MaskBelow_Button, True, False, 0)
        #Scan_ToolBox_ImageData_HBox.pack_start(Scan_ToolBox_ImageData_LoadAll_Button, True, False, 0)

        Spec_ToolBox_VBox = Gtk.VBox(homogeneous = False, spacing = 0)
        Spec_ToolBox_VBox.set_border_width(3)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_ImageData_HBox, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_CustomROI_HBox, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_CustomROI_Table, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox1, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox2, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox3, False, False, 5)
        Spec_ToolBox_VBox.set_size_request(280, 320)

        ############################  Scan Plot 1D  ###############################
        # Setup the matplotlib canvas and toolbars for 1D scan visualization

        
        self.Plot1D_Figure = Figure()
        self.Plot1D_Axe = self.Plot1D_Figure.add_axes([0.08, 0.08, 0.87, 0.87])
        self.Plot1D_Axe.set_autoscale_on(True)
        self.Plot1D_Canvas = FigureCanvas(self.Plot1D_Figure)
        Plot1D_ScrolledWindow = Gtk.ScrolledWindow()
        Plot1D_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        Plot1D_ScrolledWindow.add(self.Plot1D_Canvas)

        self.Plot1D_Canvas_Mouse_Hover_Event = self.Plot1D_Canvas.mpl_connect('motion_notify_event', self.Plot1D_Canvas_Mouse_Hover)
        self.Plot1D_Canvas_Button_Press_Event = self.Plot1D_Canvas.mpl_connect('button_press_event', self.Plot1D_Canvas_Mouse_Hover) # so that a single right click also changes the images

        self.Plot1D_LogX_ToggleButton = Gtk.ToggleButton(label = " LogX ")
        self.Plot1D_LogX_ToggleButton.connect("toggled", self.Plot1D_Vscale_Changed)
        self.Plot1D_LogY_ToggleButton = Gtk.ToggleButton(label = " LogY ")
        self.Plot1D_LogY_ToggleButton.connect("toggled", self.Plot1D_Vscale_Changed)
        Plot1D_Upload_Button = Gtk.Button(label = "+")
        Plot1D_Upload_Button.connect("clicked", self.Upload_To_Logbook, 1)
        Plot1D_Canvas_Toolbar_VSeparator1 = Gtk.VSeparator()
        Plot1D_Canvas_Toolbar_VSeparator2 = Gtk.VSeparator()

        self.Plot1D_Value_Label = Gtk.Label()     

        Plot1D_Canvas_Toolbar_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot1D_Canvas_Toolbar_HBox.set_border_width(3)
        Plot1D_Canvas_Toolbar_HBox.pack_start(Plot1D_Canvas_Toolbar_VSeparator1, False, False, 3)
        Plot1D_Canvas_Toolbar_HBox.pack_start(self.Plot1D_LogX_ToggleButton, False, False, 0)
        Plot1D_Canvas_Toolbar_HBox.pack_start(self.Plot1D_LogY_ToggleButton, False, False, 0)
        if beamline:
            Plot1D_Canvas_Toolbar_HBox.pack_start(Plot1D_Upload_Button, False, False, 0)
        Plot1D_Canvas_Toolbar_HBox.pack_start(Plot1D_Canvas_Toolbar_VSeparator2, False, False, 3)
        Plot1D_Canvas_Toolbar_HBox.pack_end(self.Plot1D_Value_Label, False, False, 0)

        Plot1D_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Plot1D_VBox.set_border_width(3)
        Plot1D_VBox.pack_start(Plot1D_Canvas_Toolbar_HBox, False, False, 0)
        Plot1D_VBox.pack_start(Plot1D_ScrolledWindow, True, True, 0)
        
        Plot1D_VBox.set_size_request(640, 480)

        #------------------------Plot 2D-------------------------------#
        # Setup the matplotlib canvas and controls (Log, Auto, Scale) for 2D scan visualization


        self.Plot2D_Figure = Figure()
        self.Plot2D_Axe = self.Plot2D_Figure.add_axes([0, 0, 1, 1])
        self.Plot2D_Canvas = FigureCanvas(self.Plot2D_Figure)
        Plot2D_ScrolledWindow = Gtk.ScrolledWindow()
        Plot2D_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        Plot2D_ScrolledWindow.add(self.Plot2D_Canvas)

        Plot2D_ScrolledWindow_EventBox = Gtk.EventBox()
        Plot2D_ScrolledWindow_EventBox.add(Plot2D_ScrolledWindow)
        
        self.Plot2D_Canvas_Mouse_Hover_Event = self.Plot2D_Canvas.mpl_connect('motion_notify_event', self.Plot2D_Canvas_Mouse_Hover)
        self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Hover) # so that a single right click also changes the images
        #self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Pressed)
        
        self.Plot2D_Log_ToggleButton = Gtk.ToggleButton(label = "Log")
        self.Plot2D_Log_ToggleButton.set_active(True)
        self.Plot2D_Log_ToggleButton.connect("toggled", self.Plot2D_Vscale_Changed)
        
        self.Plot2D_AutoScale_ToggleButton = Gtk.ToggleButton(label = "Auto")
        self.Plot2D_AutoScale_ToggleButton.set_active(True)
        self.Plot2D_AutoScale_ToggleButton.connect("toggled", self.Plot2D_AutoScale_toggled)

        self.Plot2D_ScaleBar_ToggleButton = Gtk.ToggleButton(label = "Scale")
        self.Plot2D_ScaleBar_ToggleButton.set_active(False)
        self.Plot2D_ScaleBar_ToggleButton.connect("toggled", self.Plot2D_ScaleBar_toggled)

        Plot2D_Upload_Button = Gtk.Button(label = "+")
        Plot2D_Upload_Button.connect("clicked", self.Upload_To_Logbook, 2)

        Plot2D_VSeperator = Gtk.VSeparator()

        self.Plot2D_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 10000, step_increment = .1, page_increment = 1, page_size = 0)
        self.Plot2D_Vmin_Changed_Handler = self.Plot2D_Vmin_HScale_Adjustment.connect("value_changed", self.Plot2D_Vscale_Changed)
        self.Plot2D_Vmin_HScale = Gtk.HScale()
        self.Plot2D_Vmin_HScale.set_adjustment(self.Plot2D_Vmin_HScale_Adjustment)
        self.Plot2D_Vmin_HScale.set_size_request(150, 20)
        self.Plot2D_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Plot2D_Vmin_HScale.set_digits(1)

        self.Plot2D_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 10000, lower = 10000, upper = 20000, step_increment = .1, page_increment = 1, page_size = 0)
        self.Plot2D_Vmax_Changed_Handler = self.Plot2D_Vmax_HScale_Adjustment.connect("value_changed", self.Plot2D_Vscale_Changed)
        self.Plot2D_Vmax_HScale = Gtk.HScale()
        self.Plot2D_Vmax_HScale.set_adjustment(self.Plot2D_Vmax_HScale_Adjustment)
        self.Plot2D_Vmax_HScale.set_size_request(150, 20)
        self.Plot2D_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Plot2D_Vmax_HScale.set_digits(1)
        
        self.Plot2D_Vmin_HScale.set_sensitive(False)
        self.Plot2D_Vmax_HScale.set_sensitive(False)

        self.Plot2D_P0_Label = Gtk.Label()
        self.Plot2D_P1_Label = Gtk.Label()
        
        Plot2DToolbar_HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot2DToolbar_HBox1.set_border_width(3)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Log_ToggleButton, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_AutoScale_ToggleButton, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_ScaleBar_ToggleButton, False, False, 0)
        if beamline:
           Plot2DToolbar_HBox1.pack_start(Plot2D_Upload_Button, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(Plot2D_VSeperator, False, False, 3)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Vmin_HScale, True, True, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Vmax_HScale, True, True, 0)

        Plot2DToolbar_HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot2DToolbar_HBox2.set_border_width(3)
        Plot2DToolbar_HBox2.pack_end(self.Plot2D_P0_Label, False, False, 0)
        Plot2DToolbar_HBox2.pack_start(self.Plot2D_P1_Label, False, False, 0)

        # Coarse motor position label
        self.Plot2D_CoarseMotor_Label = Gtk.Label()
        self.Plot2D_CoarseMotor_Label.set_markup(
            '<span font_desc="monospace 9" foreground="#888888">SAMX: ---  SAMY: ---  SAMZ: ---</span>')
        Plot2DToolbar_HBox3 = Gtk.HBox(homogeneous=False, spacing=3)
        Plot2DToolbar_HBox3.set_border_width(1)
        Plot2DToolbar_HBox3.pack_start(self.Plot2D_CoarseMotor_Label, False, False, 4)
 
        Plot2D_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Plot2D_VBox.set_border_width(3)
        Plot2D_VBox.pack_start(Plot2DToolbar_HBox1, False, False, 0)
        Plot2D_VBox.pack_start(Plot2DToolbar_HBox2, False, False, 0)
        Plot2D_VBox.pack_start(Plot2DToolbar_HBox3, False, False, 0)
        Plot2D_VBox.pack_start(Plot2D_ScrolledWindow_EventBox, True, True, 0)

        self.PlotSpare_Figure = Figure()
        self.PlotSpare_Axe = self.PlotSpare_Figure.add_axes([0, 0.05, 1, .9])
        self.PlotSpare_Canvas = FigureCanvas(self.PlotSpare_Figure)

        self.Plot_Notebook = Gtk.Notebook()
        self.Plot_Notebook.set_show_tabs(False)
        self.Plot_Notebook.append_page(Plot1D_VBox)
        self.Plot_Notebook.append_page(Plot2D_VBox)
        self.Plot_Notebook.append_page(self.PlotSpare_Canvas)

        #------------------------Detector Image-------------------------------#
        # Setup the main detector image display canvas and its mouse interaction events

        self.Image_Figure = Figure()
        self.Image_Axe = self.Image_Figure.add_axes([0, 0, 1, 1])
        self.Image_Axe.set_axis_off()
        #self.Image_Axe.xaxis.set_ticklabels([])
        #self.Image_Axe.yaxis.set_ticklabels([])
        self.Image_Canvas = FigureCanvas(self.Image_Figure)
        
        self.Image_Canvas_Figure_Enter_Event = self.Image_Canvas.mpl_connect('figure_enter_event', self.Image_Canvas_Figure_Entered)
        self.Image_Canvas_Button_Press_Event = self.Image_Canvas.mpl_connect('button_press_event', self.Image_Canvas_Button_Pressed)
        self.Image_Canvas_Button_Scroll_Event = self.Image_Canvas.mpl_connect('scroll_event', self.Image_Canvas_Button_Scrolled)
        self.Image_Canvas_Mouse_Hover_Event = self.Image_Canvas.mpl_connect('motion_notify_event', self.Image_Canvas_Mouse_Hover)

        self.Image_ScrolledWindow = Gtk.ScrolledWindow()
        self.Image_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.Image_ScrolledWindow.add(self.Image_Canvas)
        self.Image_ScrolledWindow_EventBox = Gtk.EventBox()
        self.Image_ScrolledWindow_EventBox.add(self.Image_ScrolledWindow)

        self.Image_Value_Label = Gtk.Label()

        self.Image_Plot_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 100000, step_increment = 1, page_increment = 1, page_size = 0)
        self.Image_Plot_HScale_Adjustment.connect("value_changed", self.Image_Plot_HScale_Changed)
        self.Image_Plot_HScale = Gtk.HScale()
        self.Image_Plot_HScale.set_adjustment(self.Image_Plot_HScale_Adjustment)
        self.Image_Plot_HScale.set_size_request(350, 20)
        self.Image_Plot_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Plot_HScale.set_digits(0)

        self.Image_ZoomOut_Button = Gtk.Button(label = "Full")
        self.Image_ZoomOut_Button.connect("clicked", self.Image_ZoomOut)
        self.Image_ZoomOut_Button.set_sensitive(False)
        self.Image_AutoScale_ToggleButton = Gtk.ToggleButton(label = "Auto")
        self.Image_AutoScale_ToggleButton.set_active(True)
        self.Image_AutoScale_ToggleButton.connect("toggled", self.Image_AutoScale_toggled)
        Image_Upload_Button = Gtk.Button(label = "+")
        Image_Upload_Button.connect("clicked", self.Upload_To_Logbook, 0)

        self.Image_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 10000, step_increment = 1, page_increment = 10, page_size = 0)
        self.Image_Vmin_Changed_Handler = self.Image_Vmin_HScale_Adjustment.connect("value_changed", self.Image_Vscale_Changed)
        self.Image_Vmin_HScale = Gtk.Scale(orientation = Gtk.Orientation.HORIZONTAL)
        self.Image_Vmin_HScale.set_adjustment(self.Image_Vmin_HScale_Adjustment)
        self.Image_Vmin_HScale.set_size_request(20, 20)
        self.Image_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Vmin_HScale.set_digits(0)
        Image_Vmin_Label = Gtk.Label()
        Image_Vmin_Label.set_text('Vmin:')

        self.Image_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 10000, lower = 10000, upper = 20000, step_increment = 1, page_increment = 10, page_size = 0)
        self.Image_Vmax_Changed_Handler = self.Image_Vmax_HScale_Adjustment.connect("value_changed", self.Image_Vscale_Changed)
        self.Image_Vmax_HScale = Gtk.HScale()
        self.Image_Vmax_HScale.set_adjustment(self.Image_Vmax_HScale_Adjustment)
        self.Image_Vmax_HScale.set_size_request(20, 20)
        self.Image_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Vmax_HScale.set_digits(0)
        Image_Vmax_Label = Gtk.Label()
        Image_Vmax_Label.set_text('Vmax:')

        self.Image_Log_ToggleButton = Gtk.ToggleButton(label = "Log")
        self.Image_Log_ToggleButton.set_active(True)
        self.Image_Log_ToggleButton.connect("toggled", self.Image_Vscale_Changed)

        self.Image_Vmin_HScale.set_sensitive(False)
        self.Image_Vmax_HScale.set_sensitive(False)
        
        # Setup toolbars for detector image manipulation (Scaling, Zoom, Sliders)
        Image_Toolbar_HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)

        Image_Toolbar_HBox1.set_border_width(3)
        Image_Toolbar_HBox1.pack_start(Image_Vmin_Label, False, False, 3)
        Image_Toolbar_HBox1.pack_start(self.Image_Vmin_HScale, True, True, 3)
        Image_Toolbar_HBox1.pack_start(Image_Vmax_Label, False, False, 3)
        Image_Toolbar_HBox1.pack_start(self.Image_Vmax_HScale, True, True, 3)

        Image_Toolbar_HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        Image_Toolbar_HBox2.set_border_width(3)
        Image_Toolbar_HBox2.pack_start(self.Image_ZoomOut_Button, False, False, 3)
        Image_Toolbar_HBox2.pack_start(self.Image_Log_ToggleButton, False, False, 3)
        Image_Toolbar_HBox2.pack_start(self.Image_AutoScale_ToggleButton, False, False, 3)
        if beamline:
            Image_Toolbar_HBox2.pack_start(Image_Upload_Button, False, False, 0)
        Image_Toolbar_HBox2.pack_start(self.Image_Plot_HScale, False, False, 0)
        Image_Toolbar_HBox2.pack_end(self.Image_Value_Label, False, False, 0)

        Image_Toolbar_VBox = Gtk.VBox()
        Image_Toolbar_VBox.set_border_width(3)
        Image_Toolbar_VBox.pack_start(Image_Toolbar_HBox1, False, False, 0)
        Image_Toolbar_VBox.pack_start(Image_Toolbar_HBox2, False, False, 0)

        Image_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Image_VBox.set_border_width(3)
        Image_VBox.pack_start(self.Image_ScrolledWindow_EventBox, True, True, 0)
        Image_VBox.pack_start(Image_Toolbar_VBox, False, False, 0)

        Image_VBox.set_size_request(740, 800)

        #------------------------XRF Plot-------------------------------#
        # Setup specialized dual-axes display for X-ray Fluorescence (XRF) spectra

        self.XRF_Figure = Figure()
        self.XRF1_Axe = self.XRF_Figure.add_axes([0.08, 0.55, 0.9, 0.43])
        self.XRF2_Axe = self.XRF_Figure.add_axes([0.08, 0.05, 0.9, 0.43])
        self.XRF_Canvas = FigureCanvas(self.XRF_Figure)
        
        #self.XRF_Canvas_Figure_Enter_Event = self.XRF_Canvas.mpl_connect('figure_enter_event', self.XRF_Canvas_Figure_Entered)
        #self.XRF_Canvas_Button_Press_Event = self.XRF_Canvas.mpl_connect('button_press_event', self.XRF_Canvas_Button_Pressed)
        #self.XRF_Canvas_Button_Scroll_Event = self.XRF_Canvas.mpl_connect('scroll_event', self.XRF_Canvas_Button_Scrolled)
        #self.XRF_Canvas_Mouse_Hover_Event = self.XRF_Canvas.mpl_connect('motion_notify_event', self.XRF_Canvas_Mouse_Hover)

        self.XRF_ScrolledWindow = Gtk.ScrolledWindow()
        self.XRF_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.XRF_ScrolledWindow.add(self.XRF_Canvas)
        self.XRF_ScrolledWindow_EventBox = Gtk.EventBox()
        self.XRF_ScrolledWindow_EventBox.add(self.XRF_ScrolledWindow)

        self.XRF_Plot_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 100000, step_increment = 1, page_increment = 1, page_size = 0)
        self.XRF_Plot_HScale_Adjustment.connect("value_changed", self.XRF_Plot_HScale_Changed)
        self.XRF_Plot_HScale = Gtk.HScale()
        self.XRF_Plot_HScale.set_adjustment(self.XRF_Plot_HScale_Adjustment)
        self.XRF_Plot_HScale.set_size_request(350, 20)
        self.XRF_Plot_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.XRF_Plot_HScale.set_digits(0)

        XRF_Sum_Button = Gtk.Button(label = "sum")
        XRF_Sum_Button.connect("clicked", self.XRF_Sum)

        XRF_RoI1_Button = Gtk.Button(label = "roi1")
        XRF_RoI1_Button.connect("clicked", self.XRF_RoI, 1)

        XRF_RoI2_Button = Gtk.Button(label = "roi2")
        XRF_RoI2_Button.connect("clicked", self.XRF_RoI, 2)

        XRF_Upload_Button = Gtk.Button(label = "+")
        XRF_Upload_Button.connect("clicked", self.Upload_To_Logbook, 3)

        self.XRF1_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 2047, step_increment = 1, page_increment = 10, page_size = 0)
        self.XRF1_Vmin_Changed_Handler = self.XRF1_Vmin_HScale_Adjustment.connect("value_changed", self.XRF1_Vscale_Changed)
        self.XRF1_Vmin_HScale = Gtk.Scale(orientation = Gtk.Orientation.HORIZONTAL)
        self.XRF1_Vmin_HScale.set_adjustment(self.XRF1_Vmin_HScale_Adjustment)
        self.XRF1_Vmin_HScale.set_size_request(20, 20)
        self.XRF1_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.XRF1_Vmin_HScale.set_digits(0)
        XRF1_Vmin_Label = Gtk.Label()
        XRF1_Vmin_Label.set_text('Vmin:')

        XRF1_Spec0_RadioButton = Gtk.RadioButton(group=None, label="1-7")
        XRF1_Spec0_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 0)
        XRF1_Spec1_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="1")
        XRF1_Spec1_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 1)
        XRF1_Spec2_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="2")
        XRF1_Spec2_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 2)
        XRF1_Spec3_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="3")
        XRF1_Spec3_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 3)
        XRF1_Spec4_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="4")
        XRF1_Spec4_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 4)
        XRF1_Spec5_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="5")
        XRF1_Spec5_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 5)
        XRF1_Spec6_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="6")
        XRF1_Spec6_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 6)
        XRF1_Spec7_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="7")
        XRF1_Spec7_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 7)
        XRF1_Spec8_RadioButton = Gtk.RadioButton(group=XRF1_Spec0_RadioButton, label="8")
        XRF1_Spec8_RadioButton.connect("toggled", self.XRF1_Spec_Changed, 8)
        self.xrf1_spec = 0

        self.XRF1_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 2047, lower = 0, upper = 2047, step_increment = 1, page_increment = 10, page_size = 0)
        self.XRF1_Vmax_Changed_Handler = self.XRF1_Vmax_HScale_Adjustment.connect("value_changed", self.XRF1_Vscale_Changed)
        self.XRF1_Vmax_HScale = Gtk.HScale()
        self.XRF1_Vmax_HScale.set_adjustment(self.XRF1_Vmax_HScale_Adjustment)
        self.XRF1_Vmax_HScale.set_size_request(20, 20)
        self.XRF1_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.XRF1_Vmax_HScale.set_digits(0)
        XRF1_Vmax_Label = Gtk.Label()
        XRF1_Vmax_Label.set_text('Vmax:')

        self.XRF2_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 2047, step_increment = 1, page_increment = 10, page_size = 0)
        self.XRF2_Vmin_Changed_Handler = self.XRF2_Vmin_HScale_Adjustment.connect("value_changed", self.XRF2_Vscale_Changed)
        self.XRF2_Vmin_HScale = Gtk.Scale(orientation = Gtk.Orientation.HORIZONTAL)
        self.XRF2_Vmin_HScale.set_adjustment(self.XRF2_Vmin_HScale_Adjustment)
        self.XRF2_Vmin_HScale.set_size_request(20, 20)
        self.XRF2_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.XRF2_Vmin_HScale.set_digits(0)
        XRF2_Vmin_Label = Gtk.Label()
        XRF2_Vmin_Label.set_text('Vmin:')

        XRF2_Spec0_RadioButton = Gtk.RadioButton(group=None, label="1-7")
        XRF2_Spec0_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 0)
        XRF2_Spec1_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="1")
        XRF2_Spec1_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 1)
        XRF2_Spec2_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="2")
        XRF2_Spec2_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 2)
        XRF2_Spec3_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="3")
        XRF2_Spec3_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 3)
        XRF2_Spec4_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="4")
        XRF2_Spec4_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 4)
        XRF2_Spec5_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="5")
        XRF2_Spec5_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 5)
        XRF2_Spec6_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="6")
        XRF2_Spec6_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 6)
        XRF2_Spec7_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="7")
        XRF2_Spec7_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 7)
        XRF2_Spec8_RadioButton = Gtk.RadioButton(group=XRF2_Spec0_RadioButton, label="8")
        XRF2_Spec8_RadioButton.connect("toggled", self.XRF2_Spec_Changed, 8)
        self.xrf2_spec = 0

        self.XRF2_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 2047, lower = 0, upper = 2047, step_increment = 1, page_increment = 10, page_size = 0)
        self.XRF2_Vmax_Changed_Handler = self.XRF2_Vmax_HScale_Adjustment.connect("value_changed", self.XRF2_Vscale_Changed)
        self.XRF2_Vmax_HScale = Gtk.HScale()
        self.XRF2_Vmax_HScale.set_adjustment(self.XRF2_Vmax_HScale_Adjustment)
        self.XRF2_Vmax_HScale.set_size_request(20, 20)
        self.XRF2_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.XRF2_Vmax_HScale.set_digits(0)
        XRF2_Vmax_Label = Gtk.Label()
        XRF2_Vmax_Label.set_text('Vmax:')

        XRF_Toolbar_HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        XRF_Toolbar_HBox1.set_border_width(3)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec0_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec1_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec2_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec3_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec4_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec5_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec6_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec7_RadioButton, False, False, 0)
        XRF_Toolbar_HBox1.pack_start(XRF1_Spec8_RadioButton, False, False, 0)
        #XRF_Toolbar_HBox1.pack_start(XRF1_Vmin_Label, False, False, 3)
        XRF_Toolbar_HBox1.pack_start(self.XRF1_Vmin_HScale, True, True, 3)
        #XRF_Toolbar_HBox1.pack_start(XRF1_Vmax_Label, False, False, 3)
        XRF_Toolbar_HBox1.pack_start(self.XRF1_Vmax_HScale, True, True, 3)

        XRF_Toolbar_HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        XRF_Toolbar_HBox2.set_border_width(3)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec0_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec1_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec2_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec3_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec4_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec5_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec6_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec7_RadioButton, False, False, 0)
        XRF_Toolbar_HBox2.pack_start(XRF2_Spec8_RadioButton, False, False, 0)
        #XRF_Toolbar_HBox2.pack_start(XRF2_Vmin_Label, False, False, 3)
        XRF_Toolbar_HBox2.pack_start(self.XRF2_Vmin_HScale, True, True, 3)
        #XRF_Toolbar_HBox2.pack_start(XRF2_Vmax_Label, False, False, 3)
        XRF_Toolbar_HBox2.pack_start(self.XRF2_Vmax_HScale, True, True, 3)

        XRF_Toolbar_HBox3 = Gtk.HBox(homogeneous = False, spacing = 3)
        XRF_Toolbar_HBox3.set_border_width(3)
        if beamline:
            XRF_Toolbar_HBox3.pack_start(XRF_Upload_Button, False, False, 0)
        XRF_Toolbar_HBox3.pack_start(XRF_Sum_Button, False, False, 0)
        XRF_Toolbar_HBox3.pack_start(XRF_RoI1_Button, False, False, 0)
        XRF_Toolbar_HBox3.pack_start(XRF_RoI2_Button, False, False, 0)
        XRF_Toolbar_HBox3.pack_end(self.XRF_Plot_HScale, False, False, 0)
        

        XRF_Toolbar_VBox = Gtk.VBox()
        XRF_Toolbar_VBox.set_border_width(3)
        XRF_Toolbar_VBox.pack_start(XRF_Toolbar_HBox1, False, False, 0)
        XRF_Toolbar_VBox.pack_start(XRF_Toolbar_HBox2, False, False, 0)
        XRF_Toolbar_VBox.pack_start(XRF_Toolbar_HBox3, False, False, 0)

        XRF_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        XRF_VBox.set_border_width(3)
        XRF_VBox.pack_start(self.XRF_ScrolledWindow_EventBox, True, True, 0)
        XRF_VBox.pack_start(XRF_Toolbar_VBox, False, False, 0)

        XRF_VBox.set_size_request(740, 800)
        ###

        # Combine Image and XRF displays into a Notebook for easy switching
        self.Image_Notebook = Gtk.Notebook()

        self.Image_Notebook.set_show_tabs(False)
        self.Image_Notebook.append_page(Image_VBox)
        self.Image_Notebook.append_page(XRF_VBox)

        # Setup the interactive python console entry
        self.Console_Entry = Gtk.Entry()

        self.Console_Entry.connect('activate', self.Console_Command_Sent)
        self.Console_Entry.connect('key-press-event', self.Console_KeyPressed)
        Console_EntryCompletion = Gtk.EntryCompletion()
        Console_EntryCompletion.set_model(Gtk.ListStore(str))
        Console_EntryCompletion.set_text_column(0)
        self.Console_Entry.set_completion(Console_EntryCompletion)
        Console_EntryCompletion.set_inline_completion(True)
        Console_EntryCompletion.set_popup_completion(False)

        # Build the application menu system (File, Settings, Add-ons, Help)
        AccelGrp = Gtk.AccelGroup.new() 
        FileMenuItem = Gtk.MenuItem(label = "File")
        FileMenu = Gtk.Menu()
        FileMenuItem.set_submenu(FileMenu)
        FileOpenMenuItem = Gtk.MenuItem(label = "Open")
        FileOpenMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_F2, 0, Gtk.AccelFlags.VISIBLE)
        FileOpenMenuItem.connect("activate", self.FileDialog_Construction, 0)
        FileRefreshMenuItem = Gtk.MenuItem(label = "Refresh")
        FileRefreshMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_F5, 0, Gtk.AccelFlags.VISIBLE)
        FileRefreshMenuItem.connect("activate", self.Folder_Refresh)
        FileQuitMenuItem = Gtk.MenuItem(label = "Quit")
        FileQuitMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_Q, Gdk.ModifierType.CONTROL_MASK, Gtk.AccelFlags.VISIBLE)
        FileQuitMenuItem.connect("activate", self.MainWindow_Destroy)
        FileMenu.add(FileOpenMenuItem)
        FileMenu.add(FileRefreshMenuItem)
        FileMenu.add(FileQuitMenuItem)
        
        SettingsMenuItem = Gtk.MenuItem(label = "Settings")
        SettingsMenu = Gtk.Menu()
        SettingsMenuItem.set_submenu(SettingsMenu)
        SettingsDetectorMenuItem = Gtk.MenuItem(label = "Detector")
        SettingsDetectorMenu = Gtk.Menu()
        SettingsDetectorMenuItem.set_submenu(SettingsDetectorMenu)
        SettingsDetector0MenuItem = Gtk.RadioMenuItem(label="1028 x 1062")
        SettingsDetector0MenuItem.set_active(True)
        SettingsDetector0MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetector1MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsDetector0MenuItem)
        SettingsDetector1MenuItem.set_label("516 x 516")
        SettingsDetector1MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetector2MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsDetector0MenuItem)
        SettingsDetector2MenuItem.set_label("515 x 515")
        SettingsDetector2MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetector3MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsDetector0MenuItem)
        SettingsDetector3MenuItem.set_label("256 x 256")
        SettingsDetector3MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetector4MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsDetector0MenuItem)
        SettingsDetector4MenuItem.set_label("487 x 195")
        SettingsDetector4MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetector5MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsDetector0MenuItem)
        SettingsDetector5MenuItem.set_label("1360 x 1024")
        SettingsDetector5MenuItem.connect("activate", self.Detector_Changed)
        SettingsDetectorMenu.add(SettingsDetector0MenuItem)
        SettingsDetectorMenu.add(SettingsDetector1MenuItem)
        SettingsDetectorMenu.add(SettingsDetector2MenuItem)
        SettingsDetectorMenu.add(SettingsDetector3MenuItem)
        SettingsDetectorMenu.add(SettingsDetector4MenuItem)
        SettingsDetectorMenu.add(SettingsDetector5MenuItem)
        SettingsColormapMenuItem = Gtk.MenuItem(label = "Colormap")
        SettingsColormapMenu = Gtk.Menu()
        SettingsColormapMenuItem.set_submenu(SettingsColormapMenu)
        SettingsColormap0MenuItem = Gtk.RadioMenuItem(label="jet")
        SettingsColormap0MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap1MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap1MenuItem.set_label("binary")
        SettingsColormap1MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap2MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap2MenuItem.set_label("gray")
        SettingsColormap2MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap3MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap3MenuItem.set_label("hot")
        SettingsColormap3MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap4MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap4MenuItem.set_label("coolwarm")
        SettingsColormap4MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap5MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap5MenuItem.set_label("viridis")
        SettingsColormap5MenuItem.set_active(True)
        SettingsColormap5MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap6MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap6MenuItem.set_label("plasma")
        SettingsColormap6MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap7MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap7MenuItem.set_label("inferno")
        SettingsColormap7MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap8MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap8MenuItem.set_label("magma")
        SettingsColormap8MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormap9MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingsColormap0MenuItem)
        SettingsColormap9MenuItem.set_label("cividis")
        SettingsColormap9MenuItem.connect("activate", self.Change_Colormap)
        SettingsColormapMenu.add(SettingsColormap0MenuItem)
        SettingsColormapMenu.add(SettingsColormap1MenuItem)
        SettingsColormapMenu.add(SettingsColormap2MenuItem)
        SettingsColormapMenu.add(SettingsColormap4MenuItem)
        SettingsColormapMenu.add(SettingsColormap5MenuItem)
        SettingsColormapMenu.add(SettingsColormap6MenuItem)
        SettingsColormapMenu.add(SettingsColormap7MenuItem)
        SettingsColormapMenu.add(SettingsColormap8MenuItem)
        SettingsColormapMenu.add(SettingsColormap9MenuItem)
        SettingsSparseItem = Gtk.CheckMenuItem(label = "Sparse Matrix")
        SettingsSparseItem.set_active(True)
        SettingsSparseItem.connect("activate", self.SparseToggled)
        SettingsShowAngleItem = Gtk.CheckMenuItem(label = "Show Angle")
        SettingsShowAngleItem.set_active(False)
        SettingsShowAngleItem.connect("activate", self.ShowAngleToggled)
        SettingsDirtyFixItem = Gtk.CheckMenuItem(label = "Dirty Fix")
        SettingsDirtyFixItem.set_active(True)
        SettingsDirtyFixItem.connect("activate", self.DirtyFixToggled)
        SettingsPumpProbeItem = Gtk.CheckMenuItem(label = "Pump Probe")
        SettingsPumpProbeItem.set_active(False)
        SettingsPumpProbeItem.connect("activate", self.PumpProbeToggled)
        SettingsXRFModeItem = Gtk.CheckMenuItem(label = "XRF Mode")
        SettingsXRFModeItem.set_active(False)
        SettingsXRFModeItem.connect("activate", self.XRFModeToggled)

        AddonsMenuItem = Gtk.MenuItem(label = "Add-ons")
        AddonsMenu = Gtk.Menu()
        AddonsMenuItem.set_submenu(AddonsMenu)
        AddonMenuItem = Gtk.MenuItem(label = "User Shortcuts")
        AddonMenuItem.connect("activate", self.Addon_Shortcuts)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Show Metadata")
        AddonMenuItem.connect("activate", self.Addon_ShowMetadata)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Run PtychoLib")
        AddonMenuItem.connect("activate", self.Addon_PtychoLib)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Probe Library")
        AddonMenuItem.connect("activate", self.Addon_ProbeLibrary)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Powder Helper")
        AddonMenuItem.connect("activate", self.Addon_PowderHelper)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Shift Correction")
        AddonMenuItem.connect("activate", self.Addon_ShiftCorrection)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Analyze 5D")
        AddonMenuItem.connect("activate", self.Addon_Analyze5D)
        AddonsMenu.add(AddonMenuItem)
            
        HelpMenuItem = Gtk.MenuItem(label = "Help")
        HelpMenu = Gtk.Menu()
        HelpMenuItem.set_submenu(HelpMenu)
        HelpAboutMenuItem = Gtk.MenuItem(label = "About")
        HelpAboutMenuItem.connect("activate", self.AboutThisProgram)
        HelpMenu.add(HelpAboutMenuItem)
        
        SettingsMenu.add(SettingsDetectorMenuItem)
        SettingsMenu.add(SettingsColormapMenuItem)
        SettingsMenu.add(SettingsSparseItem)
        SettingsMenu.add(SettingsShowAngleItem)
        SettingsMenu.add(SettingsDirtyFixItem)
        SettingsMenu.add(SettingsPumpProbeItem)
        SettingsMenu.add(SettingsXRFModeItem)
        
        MenuBar = Gtk.MenuBar()
        MenuBar.add(FileMenuItem)
        MenuBar.add(SettingsMenuItem)
        MenuBar.add(AddonsMenuItem)
        MenuBar.add(HelpMenuItem)
        
        # Finalize the main window layout and assemble all components
        self.Main_Window = Gtk.Window()

        self.Main_Window.add_accel_group(AccelGrp)

        HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox1.pack_start(MDA_File_ScrolledWindow, False, False, 0)
        HBox1.pack_start(Spec_ToolBox_VBox, False, False, 0)
        VBox1 = Gtk.VBox(homogeneous = False, spacing = 3)
        VBox1.pack_start(HBox1, False, False, 0)
        VBox1.pack_start(self.Plot_Notebook, False, False, 0)
        HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox2.pack_start(VBox1, False, False, 0)
        HBox2.pack_start(self.Image_Notebook, False, False, 0)

        # LV1
        Main_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Main_VBox.set_border_width(3)
        Main_VBox.pack_start(MenuBar, False, False, 0)
        Main_VBox.pack_start(HBox2, True, True, 0)
        Main_VBox.pack_start(self.Console_Entry, False, False, 0)

        # Main Window
        self.Main_Window.connect("destroy", self.MainWindow_Destroy)
        self.Main_Window.set_title("Diffraction User Data Explorer")
        self.Main_Window.add(Main_VBox)
        self.Main_Window.set_resizable(True)
        self.Main_Window.show_all()

        self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.TCROSS))
        self.pilatus_enabled = False
        self.eiger_enabled = True
        self.sparse_enabled= True
        self.show_angle = False
        self.dirty_fix = True
        self.pump_probe = False
        self.xrf_mode = False

#------------------------Main------------------------------#          

def main():
    # Main application loop

    
    Gtk.main()
    return 0

if __name__ == "__main__":
    # Entry point of the script


    MyMainWindow()
    main()
