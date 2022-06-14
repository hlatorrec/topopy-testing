from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFile, QgsProcessingParameterVectorDestination, QgsProcessingParameterBoolean)
from qgis import processing

from topopy import Network

class NetworkVector(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_NETW = 'INPUT_NETW'
    OUTPUT_VECTOR = 'OUTPUT_VECTOR'
    CON = 'CON'
    

    def tr(self, string):
        return QCoreApplication.translate('networkvector', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'networkvector'

    def displayName(self):
        return self.tr('Network to Vector')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        return self.tr("Load a Network object and save as a vector layer")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_NETW, self.tr('Input Network object')))
        self.addParameter(QgsProcessingParameterVectorDestination(self.OUTPUT_VECTOR, self.tr('Output vector file')))
        self.addParameter(QgsProcessingParameterBoolean(self.CON, self.tr('Split when order changes (select) or split each confluence (leave unmarked)'), False))

    def processAlgorithm(self, parameters, context, feedback):
        netw_file = self.parameterAsFile(parameters, self.INPUT_NETW, context)
        vector_file = self.parameterAsOutputLayer(parameters, self.OUTPUT_VECTOR, context)
        con = self.parameterAsBool(parameters, self.CON, context)

        netw = Network(netw_file)
        netw.export_to_shp(vector_file, con)
        results = {self.OUTPUT_VECTOR : vector_file}
        return results
