#!/usr/local/bin python3
"""Add-on for browsing end station metadata from Eiger h5 files.

Displays scalar metadata values from the beamline's end stations
(26-ID-A, 26-ID-B, 26-ID-C), Storage Ring, and other instrument
groups found under /entry/instrument/ in the h5 file.
"""

from gi.repository import Gtk
import h5py
import numpy as np

# Preferred end station to auto-select when loading a new file
_PREFERRED_STATION = "26-ID-C"

# Path inside the h5 file where end station groups live
_INSTRUMENT_PATH = "/entry/instrument"


class MainWindow:

    def Open_MetadataFile(self, widget, dude):
        """Open an h5 file via file dialog and load its metadata."""

        filename = dude.FileDialog_Construction(None, 4)
        if filename is None:
            return

        # Close any previously opened standalone h5 file
        self._close_standalone_h5()

        try:
            self.h5 = h5py.File(filename, "r")
            self._owns_h5 = True
        except Exception as e:
            print("ShowMetadata: could not open file:", e)
            return

        self.update_keystore()

    def set_h5(self, h5_file):
        """Accept an h5 file reference from dude (shared, not owned)."""

        self._close_standalone_h5()
        self.h5 = h5_file
        self._owns_h5 = False
        self.update_keystore()

    def update_keystore(self):
        """Scan the h5 file for end station groups and populate the combo box."""

        if self.h5 is None:
            return

        # Discover groups under /entry/instrument/ that contain scalar datasets
        station_keys = []
        try:
            instrument = self.h5[_INSTRUMENT_PATH]
        except KeyError:
            # Fall back: walk from root looking for groups with scalar datasets
            instrument = self.h5

        for name in instrument:
            item = instrument[name]
            if isinstance(item, h5py.Group):
                # Check if this group has at least one scalar dataset
                has_scalar = False
                for sub_name in item:
                    sub_item = item.get(sub_name)
                    if isinstance(sub_item, h5py.Dataset) and sub_item.shape == ():
                        has_scalar = True
                        break
                if has_scalar:
                    station_keys.append(name)

        # Preserve current selection if possible
        prev_active = self.Metadata_Key_ComboBox.get_active()
        prev_key = None
        if prev_active >= 0:
            model = self.Metadata_Key_ComboBox.get_model()
            if prev_active < len(model):
                prev_key = model[prev_active][1]

        self.Metadata_Key_ComboBox.handler_block(self.Key_Changed_Handler)
        self.Metadata_Key_ComboBox.get_model().clear()
        for key in station_keys:
            self.Metadata_Key_ComboBox.get_model().append([key, key])
        self.Metadata_Key_ComboBox.handler_unblock(self.Key_Changed_Handler)

        # Try to restore previous selection, or prefer 26-ID-C, or use first
        target_index = 0
        for i, key in enumerate(station_keys):
            if prev_key is not None and key == prev_key:
                target_index = i
                break
            if key == _PREFERRED_STATION:
                target_index = i

        if len(station_keys) > 0:
            self.Metadata_Key_ComboBox.set_active(target_index)

    def update(self, widget):
        """Display all scalar metadata for the selected end station."""

        self.Metadata_ListStore.clear()

        active = widget.get_active()
        if active < 0:
            return
        model = widget.get_model()
        if active >= len(model):
            return

        station_name = model[active][1]
        if self.h5 is None:
            return

        try:
            group = self.h5[_INSTRUMENT_PATH + "/" + station_name]
        except KeyError:
            try:
                group = self.h5[station_name]
            except KeyError:
                return

        for item_name in sorted(group.keys()):
            item = group.get(item_name)
            if isinstance(item, h5py.Dataset) and item.shape == ():
                raw = item[()]
                # Format the value depending on its type
                if isinstance(raw, (bytes, np.bytes_)):
                    val_str = raw.decode("utf-8", errors="replace")
                elif isinstance(raw, (float, np.floating)):
                    val_str = "{0:.6g}".format(float(raw))
                elif isinstance(raw, (int, np.integer)):
                    val_str = str(int(raw))
                else:
                    val_str = str(raw)
                self.Metadata_ListStore.append([item_name, val_str])

    def MainWindow_Destroy(self, widget, dude):
        """Clean up when the metadata window is closed."""

        self._close_standalone_h5()
        delattr(dude, "ShowMetadata")
        self.win.destroy()

    def _close_standalone_h5(self):
        """Close the h5 file only if we opened it ourselves."""

        if self._owns_h5 and self.h5 is not None:
            try:
                self.h5.close()
            except Exception:
                pass
        self.h5 = None
        self._owns_h5 = False

    def __init__(self, dude):

        self.h5 = None
        self._owns_h5 = False

        # --- End Station combo box ---
        self.Metadata_keystore = Gtk.ListStore(str, str)
        renderer_text = Gtk.CellRendererText()    
        self.Metadata_Key_ComboBox = Gtk.ComboBox.new_with_model(self.Metadata_keystore)
        self.Metadata_Key_ComboBox.pack_start(renderer_text, True)
        self.Metadata_Key_ComboBox.add_attribute(renderer_text, "text", 0)

        # --- Key/Value table ---
        Metadata_ScrolledWindow = Gtk.ScrolledWindow()
        Metadata_ScrolledWindow.set_policy(Gtk.PolicyType.ALWAYS, Gtk.PolicyType.ALWAYS)
        Metadata_ScrolledWindow.set_overlay_scrolling(False)
        self.Metadata_ListStore = Gtk.ListStore(str, str)
        self.Metadata_TreeView_Filter = self.Metadata_ListStore.filter_new()
        self.Metadata_TreeView = Gtk.TreeView.new_with_model(self.Metadata_TreeView_Filter)
        self.Metadata_TreeView.get_selection().set_mode(Gtk.SelectionMode.NONE)

        Metadata_CellRendererText = Gtk.CellRendererText()     
        #Metadata_CellRendererText.set_property("xalign", 1)
        Metadata_TreeViewColumn = Gtk.TreeViewColumn("key", Metadata_CellRendererText, text = 0)
        self.Metadata_TreeView.append_column(Metadata_TreeViewColumn)

        Metadata_CellRendererText = Gtk.CellRendererText()     
        #Metadata_CellRendererText.set_property("xalign", 1)
        Metadata_TreeViewColumn = Gtk.TreeViewColumn("value", Metadata_CellRendererText, text = 1)
        self.Metadata_TreeView.append_column(Metadata_TreeViewColumn)

        self.Key_Changed_Handler = self.Metadata_Key_ComboBox.connect("changed", self.update)

        self.Metadata_TreeView.set_headers_visible(False)
        self.Metadata_TreeView.set_enable_search(False)
        Metadata_ScrolledWindow.add(self.Metadata_TreeView)
        Metadata_ScrolledWindow.set_size_request(300, 400)

        # --- Open button ---
        Metadata_FileOpen_Button = Gtk.Button(" Open ")
        Metadata_FileOpen_Button.connect("clicked", self.Open_MetadataFile, dude)
        
        # --- Layout ---
        Main_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Main_HBox.set_border_width(3)
        Main_HBox.pack_start(Metadata_FileOpen_Button, False, False, 0)
        Main_HBox.pack_start(self.Metadata_Key_ComboBox, True, True, 0)

        Main_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Main_VBox.set_border_width(3)
        Main_VBox.pack_start(Main_HBox, False, False, 0)
        Main_VBox.pack_start(Metadata_ScrolledWindow, False, False, 0)

        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        self.win.add(Main_VBox)
        self.win.set_title("Show Metadata")
        self.win.show_all()

        # Auto-load from dude's currently open h5 file if available
        if hasattr(dude, "h5") and dude.h5 is not None:
            try:
                dude.h5.id  # verify the file is still open
                self.h5 = dude.h5
                self._owns_h5 = False
                self.update_keystore()
            except Exception:
                pass


def main():
    
    Gtk.main()
    return 0

if __name__ == "__main__":

    MyMainWindow()
    main()
