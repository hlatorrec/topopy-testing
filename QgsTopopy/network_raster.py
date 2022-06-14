from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFile, QgsProcessingParameterRasterDestination)
from qgis import processing

from topopy import Network

class NetworkRaster(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_NETW = 'INPUT_NETW'
    OUTPUT_RASTER = 'OUTPUT_RASTER'

    def tr(self, string):
        return QCoreApplication.translate('networkraster', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'networkraster'

    def displayName(self):
        return self.tr('Network to Raster')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        return self.tr("Load a Network object and save as a raster layer")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_NETW, self.tr('Input Network object')))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_RASTER, self.tr('Output network raster'), None, False))

    def processAlgorithm(self, parameters, context, feedback):
        netw_file = self.parameterAsFile(parameters, self.INPUT_NETW, context)
        netw_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_RASTER, context)
        
        netw = Network(netw_file)
        netw_as_raster = netw.get_streams()
        netw_as_raster.save(netw_raster)
        results = {self.OUTPUT_RASTER : netw_raster}
        return results