###############################################################################

.AnnotationData <- setClass("AnnotationData",
        representation(varMetadata = "data.frame",
                data = "data.frame",
                dimLabels = "character"),
        prototype = prototype(
                varMetadata = new( "data.frame" ),
                data = new( "data.frame" ),
                dimLabels=c("rowNames", "columnNames")))
           
setClass("AnnotatedData",
        representation(assayData = "matrix",
                phenoData = "AnnotationData",
                featureData = "AnnotationData"),
        prototype = prototype(
                phenoData = .AnnotationData(
                        dimLabels=c("sampleNames", "sampleColumns")),
                featureData = .AnnotationData(
                        dimLabels=c("featureNames", "featureColumns"))
        )
)

# generics
setGeneric('colData', function(object) standardGeneric('colData'))
setGeneric('colData<-', function(object, value) standardGeneric('colData<-'))
#
setGeneric('rowData', function(object) standardGeneric('rowData'))
setGeneric('rowData<-', function(object, value) standardGeneric('rowData<-'))
#
setGeneric('mainData', function(object) standardGeneric('mainData'))
setGeneric('mainData<-', function(object, value) standardGeneric('mainData<-'))
#
setGeneric('varData', function(object) standardGeneric('varData'))
