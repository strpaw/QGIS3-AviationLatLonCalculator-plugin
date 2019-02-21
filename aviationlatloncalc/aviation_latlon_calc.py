# -*- coding: utf-8 -*-
"""
/***************************************************************************
 AviationLatLonCalc
                                 A QGIS plugin
 Plugin to calculate latitude and longitude when loacitons given in aviation lake way (azimuth, bearings, offesets against reference point)
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2019-02-16
        git sha              : $Format:%H$
        copyright            : (C) 2019 by Paweł Strzelewicz
        email                : aviationgisapp@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt5.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, QVariant
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QMessageBox, QWidget
from qgis.core import *

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .aviation_latlon_calc_dialog import AviationLatLonCalcDialog
import os.path
from . import aviation_gis_core as agc

w = QWidget()


class AviationLatLonCalc:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        self.ref_point = None
        self.mlyr_name = ''
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'AviationLatLonCalc_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = AviationLatLonCalcDialog()
        self.dlg.comboBoxInputData.currentIndexChanged.connect(self.set_calculated_point_input)
        self.dlg.pushButtonCalculate.clicked.connect(self.calculate)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&AviationLatLonCalculator')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'AviationLatLonCalc')
        self.toolbar.setObjectName(u'AviationLatLonCalc')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('AviationLatLonCalc', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/aviation_latlon_calc/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'AviationLatLonCalculator'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&AviationLatLonCalculator'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    @staticmethod
    def create_mem_lyr(lyr_name):
        """ Create temporary 'memory' layer to store results.
        :param lyr_name: string, layer name
        """
        mlyr = QgsVectorLayer('Point?crs=epsg:4326', lyr_name, 'memory')
        prov = mlyr.dataProvider()
        mlyr.startEditing()
        prov.addAttributes([QgsField("CP_ID", QVariant.String),  # Origin point ID
                            QgsField("CP_LATDD", QVariant.String),  # Calculated point latitude (decimal degrees)
                            QgsField("CP_LONDD", QVariant.String),  # Calculated point longitude (decimal degrees)
                            QgsField("CP_DEF", QVariant.String)])  # Calculated point in polar coordinates
        mlyr.commitChanges()
        QgsProject.instance().addMapLayer(mlyr)

    def add_point_to_layer(self, point, attributes):
        out_lyr = self.iface.activeLayer()
        out_lyr.startEditing()
        out_prov = out_lyr.dataProvider()
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPointXY(point))
        feat.setAttributes(attributes)
        out_prov.addFeatures([feat])
        out_lyr.commitChanges()
        out_lyr.updateExtents()
        self.iface.mapCanvas().setExtent(out_lyr.extent())
        self.iface.mapCanvas().refresh()

    def set_calculated_point_input(self):
        """ Sets controls to enter data for calculated point according to chosen method/input data """
        if self.dlg.comboBoxInputData.currentIndex() == 0:  # Azimuth, distance
            self.dlg.stackedWidgetCalcPointInput.setCurrentIndex(0)
        elif self.dlg.comboBoxInputData.currentIndex() == 1:  # Azimuth, distance, offset
            self.dlg.stackedWidgetCalcPointInput.setCurrentIndex(1)
        elif self.dlg.comboBoxInputData.currentIndex() == 2:  # Local Cartesian
            self.dlg.stackedWidgetCalcPointInput.setCurrentIndex(2)

    def check_ref_point(self):
        """ Checks if reference point input data is correct
        :return: check_result: bool, True if is correct, False otherwise
        :return: check_msg: str, empty string if input is correct, string with errors otherwise
        """
        check_result = True
        check_msg = ''
        self.ref_point = agc.AviationBasePoint(self.dlg.lineEditID.text(),
                                               self.dlg.lineEditRefLat.text(),
                                               self.dlg.lineEditRefLon.text(),
                                               self.dlg.lineEditRefMagVar.text())
        if self.ref_point.is_valid is False:
            check_result = False
            check_msg += self.ref_point.err_msg

        return check_result, check_msg

    def get_ad_distance_uom(self):
        """ Gets unit of measure for azimuth, distance single point method """
        if self.dlg.comboBoxADPointDistanceUOM.currentIndex() == 0:
            return agc.UOM_M
        if self.dlg.comboBoxADPointDistanceUOM.currentIndex() == 1:
            return agc.UOM_KM
        if self.dlg.comboBoxADPointDistanceUOM.currentIndex() == 2:
            return agc.UOM_NM
        if self.dlg.comboBoxADPointDistanceUOM.currentIndex() == 3:
            return agc.UOM_FEET
        if self.dlg.comboBoxADPointDistanceUOM.currentIndex() == 4:
            return agc.UOM_SM

    def get_ado_distance_uom(self):
        """ Gets unit of measure for azimuth, distance, offset single point method """
        if self.dlg.comboBoxADOPointDistanceUOM.currentIndex() == 0:
            return agc.UOM_M
        if self.dlg.comboBoxADOPointDistanceUOM.currentIndex() == 1:
            return agc.UOM_KM
        if self.dlg.comboBoxADOPointDistanceUOM.currentIndex() == 2:
            return agc.UOM_NM
        if self.dlg.comboBoxADOPointDistanceUOM.currentIndex() == 3:
            return agc.UOM_FEET
        if self.dlg.comboBoxADOPointDistanceUOM.currentIndex() == 4:
            return agc.UOM_SM

    def get_ado_offset_side(self):
        if self.dlg.comboBoxADOPointOffsetSite.currentIndex() == 0:
            return 'LEFT'
        elif self.dlg.comboBoxADOPointOffsetSite.currentIndex() == 1:
            return 'RIGHT'

    def get_ado_offset_uom(self):
        if self.dlg.comboBoxADOOffsetPointDistanceUOM.currentIndex() == 0:
            return agc.UOM_M
        if self.dlg.comboBoxADOOffsetPointDistanceUOM.currentIndex() == 1:
            return agc.UOM_KM
        if self.dlg.comboBoxADOOffsetPointDistanceUOM.currentIndex() == 2:
            return agc.UOM_NM
        if self.dlg.comboBoxADOOffsetPointDistanceUOM.currentIndex() == 3:
            return agc.UOM_FEET
        if self.dlg.comboBoxADOOffsetPointDistanceUOM.currentIndex() == 4:
            return agc.UOM_SM

    def get_cartesian_x_uom(self):
        if self.dlg.comboBoxCartesianXUOM.currentIndex() == 0:
            return agc.UOM_M
        if self.dlg.comboBoxCartesianXUOM.currentIndex() == 1:
            return agc.UOM_KM
        if self.dlg.comboBoxCartesianXUOM.currentIndex() == 2:
            return agc.UOM_NM
        if self.dlg.comboBoxCartesianXUOM.currentIndex() == 3:
            return agc.UOM_FEET
        if self.dlg.comboBoxCartesianXUOM.currentIndex() == 4:
            return agc.UOM_SM

    def get_cartesian_y_uom(self):
        if self.dlg.comboBoxCartesianYUOM.currentIndex() == 0:
            return agc.UOM_M
        if self.dlg.comboBoxCartesianYUOM.currentIndex() == 1:
            return agc.UOM_KM
        if self.dlg.comboBoxCartesianYUOM.currentIndex() == 2:
            return agc.UOM_NM
        if self.dlg.comboBoxCartesianYUOM.currentIndex() == 3:
            return agc.UOM_FEET
        if self.dlg.comboBoxCartesianYUOM.currentIndex() == 4:
            return agc.UOM_SM

    def calculate(self):
        """ Checks if input data for calculated point is correct """
        check_result, check_msg = self.check_ref_point()
        calc_point_id = ''

        calc_point = agc.AviationCalculatedPoint(self.ref_point)

        if self.dlg.comboBoxInputData.currentIndex() == 0:  # Azimuth, Distance single point
            ad_brng = agc.Bearing(self.dlg.lineEditADPointAzimuth.text())
            ad_dist = agc.Distance(self.dlg.lineEditADPointDistance.text(),
                                   self.get_ad_distance_uom())

            if ad_brng.is_valid is False:
                check_result = False
                check_msg += ad_brng.err_msg

            if ad_dist.is_valid is False:
                check_result = False
                check_msg += ad_dist.err_msg

            if check_result is False:
                QMessageBox.critical(w, "Message", check_msg)
            else:
                calc_point_id = self.dlg.lineEditADPointID.text()
                calc_point.polar_coordinates2latlon(ad_brng, ad_dist)

        elif self.dlg.comboBoxInputData.currentIndex() == 1:
            ado_distance = agc.Distance(self.dlg.lineEditADOPointDistance.text(),
                                        self.get_ado_distance_uom())
            ado_brng = agc.Bearing(self.dlg.lineEditADOPointAzimuth.text())
            ado_offset_side = self.get_ado_offset_side()
            ado_offset_distance = agc.Distance(self.dlg.lineEditADOOffsetPointDistance.text(),
                                               self.get_ado_offset_uom(), checked_value='Offset distance')

            if ado_distance.is_valid is False:
                check_result = False
                check_msg += ado_distance.err_msg

            if ado_brng.is_valid is False:
                check_result = False
                check_msg += ado_brng.err_msg

            if ado_offset_distance.is_valid is False:
                check_result = False
                check_msg += ado_offset_distance.err_msg

            if check_result is False:
                QMessageBox.critical(w, "Message", check_msg)
            else:
                calc_point_id = self.dlg.lineEditADOPointID.text()
                calc_point.offset_coordinates2latlon(ado_brng, ado_distance, ado_offset_side, ado_offset_distance)

        elif self.dlg.comboBoxInputData.currentIndex() == 2:

            x_axis_brng = agc.Bearing(self.dlg.lineEditXAxisBearing.text(), checked_value='X axis bearing')
            x = agc.Distance(self.dlg.lineEditCartesianX.text(), self.get_cartesian_x_uom(),
                             checked_value='X coordinate', allow_negative=True)
            y = agc.Distance(self.dlg.lineEditCartesianY.text(), self.get_cartesian_y_uom(),
                             checked_value='Y coordinate', allow_negative=True)

            if x_axis_brng.is_valid is False:
                check_result = False
                check_msg += x_axis_brng.err_msg

            if x.is_valid is False:
                check_result = False
                check_msg += x.err_msg

            if y.is_valid is False:
                check_result = False
                check_msg += y.err_msg

            if check_result is False:
                QMessageBox.critical(w, "Message", check_msg)
            else:
                calc_point_id = self.dlg.lineEditLCPointID.text()
                y_axis_orient = self.dlg.comboBoxYAxisOrientation.currentText()
                calc_point.local_cartesian_coordinates2latlon(x_axis_brng, y_axis_orient, x, y)

        if check_result is True:
            cp_qgs_point = QgsPointXY(calc_point.cp_lon_dd, calc_point.cp_lat_dd)
            cp_attributes = [calc_point_id,
                             str(calc_point.cp_lat_dd),
                             str(calc_point.cp_lon_dd),
                             calc_point.cp_definition]

            layers = QgsProject.instance().layerTreeRoot().children()
            layers_list = []  # List of layers in current (opened) QGIS project
            for layer in layers:
                layers_list.append(layer.name())

            if self.mlyr_name == '':
                self.mlyr_name = 'TEST'

            if self.mlyr_name not in layers_list:
                self.create_mem_lyr(self.mlyr_name)
                self.add_point_to_layer(cp_qgs_point, cp_attributes)
            else:
                self.add_point_to_layer(cp_qgs_point, cp_attributes)

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass