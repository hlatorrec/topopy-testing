from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsFeatureSink, QgsProcessingException, QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterRasterDestination, QgsProcessingParameterField, QgsProcessingParameterField)
from qgis import processing

from topopy import Flow, extract_points
from numpy import array

class DrainageBasin(QgsProcessingAlgorithm):
    INPUT_FD = 'INPUT_FD'
    OUTPUT_DBAS = 'OUTPUT_DBAS'
    OUTLETS = 'OUTLETS'
    ID_FIELD = 'ID_FIELD'

    def tr(self, string):
        return QCoreApplication.translate('drainagebasin', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'drainagebasin'

    def displayName(self):
        return self.tr('Drainage Basin')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        texto = 'Placeholder text'
        return texto

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_FD,  self.tr("Input Flow Direction raster")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_DBAS, "Output Drainage Basin raster", None, False))
        self.addParameter(QgsProcessingParameterFeatureSource(self.OUTLETS, 'Outlets', [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField(self.ID_FIELD, 'Id Field', None, self.OUTLETS, QgsProcessingParameterField.Numeric))

    def processAlgorithm(self, parameters, context, feedback):
        fd_raster = self.parameterAsRasterLayer(parameters, self.INPUT_FD, context)
        dbas_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_DBAS, context)
        outlets = self.parameterAsSource(parameters, self.OUTLETS, context)
        
        fd = Flow(fd_raster.source()) #Load Flow Direction raster
        #Get 0.25 threshold area
        threshold = int(fd.get_ncells() * 0.0025)
        id_field = self.parameterAsString(parameters, self.ID_FIELD, context)
        points = []
        for feat in outlets.getFeatures():
            geom = feat.geometry().asPoint()
            oid = feat.attribute(id_field)
            points.append([geom.x(), geom.y(), oid])
        points = array(points)
        outlets = fd.snap_points(points, threshold, 'channel', True)
        dbas = fd.get_drainage_basins(outlets=outlets)
        dbas.save(dbas_raster) #Save drainage basins raster onto disc
        results = {self.OUTPUT_DBAS : dbas_raster}
        return results
