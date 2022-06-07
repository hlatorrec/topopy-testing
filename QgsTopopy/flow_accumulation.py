from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsFeatureSink, QgsProcessingException, QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterBoolean, QgsProcessingParameterRasterDestination)
from qgis import processing

from topopy import Flow

class FlowAccumulation(QgsProcessingAlgorithm):
    INPUT_FD = 'INPUT_FD' #Input DEM grid
    OUTPUT_FAC = 'OUTPUT_FAC' #Output flow accumulation grid
        
    def tr(self, string):
        return QCoreApplication.translate('flowaccumulation', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'flowaccumulation'

    def displayName(self):
        return self.tr('Flow Accumulation')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        texto = "Takes a DEM as inpuut and outputs a flow direction raster using the Flow.get_flow_accumulation() method. Keeps default parameters (mantains NoData values, outputs a topopy.Grid object)."
        return texto

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_FD,  self.tr("Input Flow Direction raster")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FAC, "Output Flow Accumulation raster", None, False))
        
    def processAlgorithm(self, parameters, context, feedback):
        fd_raster = self.parameterAsRasterLayer(parameters, self.INPUT_FD, context)
        fac_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_FAC, context)
        
        fd = Flow(fd_raster.source()) #Load Flow object
        fac = fd.get_flow_accumulation() #Get flow accumulation raster
        fac.save(fac_raster) #Save flow accumulation raster onto disc
        results = {self.OUTPUT_FAC : fac_raster}
        return results
