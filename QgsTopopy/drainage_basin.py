from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsFeatureSink, QgsProcessingException, QgsProcessingAlgorithm, QgsProcessingParameterFeatureSource, QgsProcessingParameterFeatureSink, QgsProcessingParameterRasterLayer, parameterAsOutputLayer, QgsProcessingParameterNumber)
from qgis import processing

from topopy import DEM, Flow

class DrainageBasin(QgsProcessingAlgorithm):
    INPUT_DEM = 'INPUT_DEM'
    OUTPUT_DBAS = 'OUTPUT_DBAS'
    OUTLETS = 'OUTLETS'
    MIN_AREA = 'MIN_AREA'

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
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM,  self.tr("Input DEM")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_DBAS, "Output Drainage Basin raster", None, False))
        self.addParameter(QgsProcessingParameterNumber(self.MIN_AREA, "Minimum area for basins (avoids very small basins)", type=QgsProcessingParameterNumber.Double))
        
    def processAlgorithm(self, parameters, context, feedback):
        dem_raster = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        dbas_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_DBAS, context)
        min_area = self.parameterAsDouble(parameters, self.MIN_AREA, context)
        
        dem = DEM(dem_raster.source()) #Load DEMProcessing
        fd = Flow(dem)
        dbas = fd.get_drainage_basins(min_area = min_area)
        dbas.save(dbas_raster) #Save drainage basins raster onto disc
        results = {self.OUTPUT_DBAS : dbas_raster}
        return results