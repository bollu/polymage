import ctypes
import _ctypes
import impipe
import constructs
import subprocess
import time
import string
import random
import os
from common import clock, StatValue, printLine

# TODO: 

# 1. Decide on all required parameters to the 'generate' and 'execute'  ( )
#    functions.
#

# TODO:
# 1. Design in a robust way to handle encompassing of various search    ( )
#    spaces.
# 2. Introduce parallelism in code generation and compilation           ( )
# 3. Make the search configurations in each space a Set before 
#    enumerating                                                        ( )
#
def generate(_pm_argDict):

    # unpack the arguments from the arg dictionary
    try:
        _pm_appName = _pm_argDict['_pm_appName']
    except KeyError:
        print 'tuner : generator : \'_pm_appName\' - \
                not an optional parameter'

    try:
        _pm_liveOuts = _pm_argDict['_pm_liveOuts']
    except KeyError:
        print 'tuner : generator : \'_pm_liveOuts\' - \
                not an optional parameter'

    try:
        _pm_paramConstraints = _pm_argDict['_pm_paramConstraints']
    except KeyError:
        _pm_paramConstraints = None

    try:
        _pm_paramEstimates = _pm_argDict['_pm_paramEstimates']
    except KeyError:
        _pm_paramEstimates = None

    try:
        _pm_tileSizeConfigs = _pm_argDict['_pm_tileSizeConfigs']
    except KeyError:
        _pm_tileSizeConfigs = None

    try:
        _pm_groupSizeConfigs = _pm_argDict['_pm_groupSizeConfigs']
    except KeyError:
        _pm_groupSizeConfigs = None

    # for now, assume the following option to be deprecated
    '''
    try:
        _pm_inlineDirectives = _pm_argDict['_pm_inlineDirectives']
    except KeyError:
        _pm_inlineDirectives = None
    '''

    try:
        _pm_dstPath = _pm_argDict['_pm_dstPath']
    except KeyError:
        _pm_dstPath = '/tmp'

    try:
        _pm_shouldDebug = _pm_argDict['_pm_shouldDebug']
    except KeyError:
        _pm_shouldDebug = False
   


    # TODO:
    # 1. write a neat printer                   ( )
    #
    def printParams(toFile):
        if toFile:
            printLine(toFile)
            print >>toFile, "App name                : \""+_pm_appName+"\""
            print >>toFile, "Code Path               : \""+_pm_dstPath+"\""
            print >>toFile, "Live-out Functions      :",
            for func in _pm_liveOuts:
                print >>toFile, "\""+func.name+"\",",
            print >>toFile, "\nParameter Constraints   :"
            for constraint in _pm_paramConstraints:
                print >>toFile, constraint
            print >>toFile, "\nParameter Estimates     :"
            for estimate in _pm_paramEstimates:
                print >>toFile, (estimate[0].name, estimate[1])
            print >>toFile, "\nTile Sizes              :"
            for tileSizes in _pm_tileSizeConfigs:
                print >>toFile, tileSizes
            print >>toFile, "\nGroup Sizes             :"
            for groupSizes in _pm_groupSizeConfigs:
                print >>toFile, groupSizes
        if _pm_shouldDebug:
            printLine()
            print "App name                : \""+_pm_appName+"\""
            print "Code Path               : \""+_pm_dstPath+"\""
            print "Live-out Functions      :",
            for func in _pm_liveOuts:
                print "\""+func.name+"\",",
            print
            print "\nParameter Constraints :"
            for constraint in _pm_paramConstraints:
                print constraint
            print "\nParameter Estimates   :"
            for estimate in _pm_paramEstimates:
                print (estimate[0].name, estimate[1])
            print "\nTile Sizes            :"
            for tileSizes in _pm_tileSizeConfigs:
                print tileSizes
            print "\nGroup Sizes           :"
            for groupSizes in _pm_groupSizeConfigs:
                print groupSizes
            '''
            print "\nfunctions to be inlined:"
            for func in _pm_inlineDirectives:
                print "\""+func.name+"\",",
            '''

    appName = _pm_appName+'_polymage_'

    def randomString():
        # for now, using a random string of length 10 
        # as subdir name
        return ''.join(random.SystemRandom().choice(
                    string.lowercase + \
                    string.uppercase + \
                    string.digits) \
                    for _ in xrange(10))
 
    # ensure that a subdirectory is created, with a name not conflicting with
    # existing ones
    dstSubDir = str(_pm_dstPath)+'/'+'Poly'+randomString()+'Mage/'
    while os.path.exists(dstSubDir):
        dstSubDir =str(_pm_dstPath)+'/'+'Poly'+randomString()+'Mage/'

    # subdirectories and files
    os.makedirs(dstSubDir)
    progPrefix = str(dstSubDir)+str(appName)
    configFileName = str(dstSubDir)+'configurations.txt'
    configFile = open(configFileName, 'w')

    # Compile String parts
    cxx='icpc'
    #cxx='g++-4.8'
    optFlags='-openmp -xhost -O3 -ipo'
    #optFlags='-fopenmp -O3'
    sharedLibFlags='-fPIC -shared'

    # Generate group sizes automatically, if none is specified by the user.
    # TODO:
    # 1. Limit the configuration space to the total number of
    #    functions in the pipeline, since any group size greater
    #    than that will lead to redundancy.                             ( )
    #
    if _pm_groupSizeConfigs == None:
        for i in range(1, 4):
            _pm_groupSizeConfigs.append(i)
        # Assuming that the number of functions does not exceed
        # this maximum, include a "group-all"
        _pm_groupSizeConfigs.append(200)

    '''
    # Generate tile sizes automatically, if none are specified by the user.
    # However, it makes no sense to iterate over tileSize space for 4D when
    # the pipeline actually has lesser dimensionality. 

    # TODO:
    # 1. This has to be further improved to use DSL compiler 
    #    information
    #    -- restrict to tileable dimensions.                            ( )
    #    -- target multiple architectures etc.                          ( )
    #    -- possible to restrict tileSize space? (based on 
    #       functions size)                                             ( )

    if _pm_tileSizeConfigs == None:
        tileBase1 = []
        for tile1 in [8, 16, 32, 64]:
            tileBase2 = tileBase1+[tile1]
            for tile2 in [8, 16, 32, 64, 128]:
                tileBase3 = tileBase2+[tile2]
                for tile3 in [8, 16, 32, 64, 128, 256]:
                    tileBase4 = tileBase3+[tile3]
                    for tile4 in [8, 16, 32, 64, 128, 256]:
                        tileSize = tileBase4+[tile4]
                        print tileSize
                        _pm_tileSizeConfigs.append(tileSize)
    '''

    printParams(configFile)

    # TODO:
    # 1. modify this whole looping, when 'thresholding' is implemented
    #    in the DSL compiler                                            ( )
    #
    _pm_config = 0

    total_t1 = clock()
    # iterate over groupSizes
    for _pm_groupSize in _pm_groupSizeConfigs:
        for _pm_tileSize in _pm_tileSizeConfigs:
            _pm_config += 1

            # Update configs file:
            printLine(configFile)
            print >>configFile, "Config     : #"+str(_pm_config)
            print >>configFile, "Tile Sizes : "+str(_pm_tileSize)
            print >>configFile, "Group Size : "+str(_pm_groupSize)

            # .cpp and .so files
            cFileName = str(progPrefix)+str(_pm_config)+'.cpp'
            soFileName = str(progPrefix)+str(_pm_config)+'.so'

            t1 = clock()

            # building the pipeline :
            _pm_pipe = impipe.buildPipeLine(
                            _pm_liveOuts,
                            paramConstraints=_pm_paramConstraints,
                            paramEstimates=_pm_paramEstimates,
                            tileSizes=_pm_tileSize,
                            groupSize=_pm_groupSize)
                            #inlineDirectives=_pm_inlineDirectives)
 
            # code generation :
            cFile = open(cFileName, 'w')
            cFile.write(_pm_pipe.generateCNaive(isExternAlloc=True, \
                                                isExternFunc=True, \
                                                areParamsVoidPtrs=True).__str__())
            cFile.close()

            t2 = clock()

            codegenTime = float(t2) - float(t1)
            print >>configFile, "Code Generation Time   : "+str(codegenTime*1000)

            # TODO:
            # 1. Currently inlineDirectives changes the whole DAG
            #    of the pipeline without cloning. Add this parameter
            #    after the compiler facilitates modification.           ( )
            #
  
            # TODO:
            # 1. Handle other CXX flags                                 ( )
            # 2. Include g++ support                                    ( )
            # 3. Should this be using just a 'make'?                    ( )
            #
            compileString = cxx+' '+ \
                            optFlags+' '+ \
                            sharedLibFlags+' '+ \
                            cFileName+' -o'+ \
                            soFileName

            t1 = clock()

            # compilation :
            try:
                subprocess.check_output(compileString, shell=True)
                pass
            finally:
                pass

            t2 = clock()

            compileTime = float(t2) - float(t1)
            print >>configFile, "Compilation Time       : "+str(compileTime*1000)

            # total time for this variant:
            print >>configFile, "Total                  : "+str((codegenTime+compileTime)*1000)

    total_t2 = clock()
    totalTime = float(total_t2) - float(total_t1)

    printLine(configFile)
    print >>configFile, "Time taken to generate all variants : "
    print >>configFile, str(totalTime*1000)

    printLine(configFile)
    configFile.close()

    return dstSubDir, _pm_config

def execute(_pm_argDict):

    # unpack the arguments from the arg dictionary
    try:
        _pm_pipeArgDict = _pm_argDict['_pm_pipeArgDict']
    except KeyError:
        print 'tuner : executer : \'_pm_pipeArgDict\' - \
                not an optional parameter'

    try:
        _pm_appName = _pm_argDict['_pm_appName']
    except KeyError:
        print 'tuner : executer : \'_pm_appName\' - \
                not an optional parameter'

    try:
        _pm_pipe = _pm_argDict['_pm_pipe']
    except KeyError:
        print 'tuner : executer : \'_pm_pipe\' - \
                not an optional parameter'

    try:
        _pm_srcPath = _pm_argDict['_pm_srcPath']
    except KeyError:
        _pm_srcPath = "/tmp"

    try:
        _pm_configsCount = _pm_argDict['_pm_configsCount']
    except KeyError:
        _pm_configsCount = 0

    try:
        _pm_omp_threads = _pm_argDict['_pm_omp_threads']
    except KeyError:
        _pm_omp_threads = 1

    try:
        _pm_nruns = _pm_argDict['_pm_nruns']
    except KeyError:
        _pm_nruns = 1

    try:
        _pm_shouldDebug = _pm_argDict['_pm_shouldDebug']
    except KeyError:
        _pm_shouldDebug = True

    def printParams(toFile=None):
        if toFile:
            printLine(toFile)
            print >>toFile, "App Name              : \""+_pm_appName+"\""
            print >>toFile, "Total Configurations  :", _pm_configsCount
            print >>toFile, "OMP_NUM_THREADS       :",_pm_omp_threads
            print >>toFile, "Number of Tuning Runs :",_pm_nruns
        if _pm_shouldDebug:
            printLine()
            print "app name              : \""+_pm_appName+"\""
            print "Total Configurations  :", _pm_configsCount
            print "OMP_NUM_THREADS       :",_pm_omp_threads
            print "Number of tuning runs :",_pm_nruns

    # Parameters
    _pm_pipeParams = _pm_pipe.getParameters()
    _pm_pipeParams.sort(key=lambda x: x.name)

    # Inputs (Images)
    _pm_pipeInputs = _pm_pipe.inputs
    _pm_pipeInputs.sort(key=lambda x: x.name)

    # Outputs
    _pm_pipeOutputs = _pm_pipe.outputs
    _pm_pipeOutputs.sort(key=lambda x: x.name)

    # construct shared library function properties

    def convertToCType(inpType, inpValue):
        if inpType == 'void':
            return ctypes.c_void(inpValue)

        if inpType == 'char':
            return ctypes.c_char(inpValue)
        if inpType == 'unsigned char':
            return ctypes.c_ubyte(inpValue)

        if inpType == 'short' or \
            inpType == 'short int':
            return ctypes.c_short(inpValue)
        if inpType == 'unsigned short' or \
            inpType == 'unsigned short int':
            return ctypes.c_ushort(inpValue)

        if inpType == 'int':
            return ctypes.c_int(inpValue)
        if inpType == 'unsigned' or \
            inpType == 'unsigned int':
            return ctypes.c_uint(inpValue)

        if inpType == 'long' or \
            inpType == 'long int':
            return ctypes.c_long(inpValue)
        if inpType == 'unsigned long' or \
            inpType == 'unsigned long int':
            return ctypes.c_ulong(inpValue)

        if inpType == 'long long' or \
            inpType == 'long long int':
            return ctypes.c_longlong(inpValue)
        if inpType == 'unsigned long long' or \
            inpType == 'unsigned long long int':
            return ctypes.c_ulonglong(inpValue)

        if inpType == 'float':
            return ctypes.c_float(inpValue)
        if inpType == 'double':
            return ctypes.c_double(inpValue)
        if inpType == 'long double':
            return ctypes.c_double(inpValue)

    libFunctionName = 'pipeline_'+_pm_pipe.name

    args = []

    for param in _pm_pipeParams:
        args += [convertToCType(param.typ().cName(), _pm_pipeArgDict[param.name])]

    for inp in _pm_pipeInputs:
        args += [ctypes.c_void_p(_pm_pipeArgDict[inp.name].ctypes.data)]

    for out in _pm_pipeOutputs:
        args += [ctypes.c_void_p(_pm_pipeArgDict[out.name].ctypes.data)]


    # set other variables
    appName = _pm_appName+'_polymage_'
    progPrefix = str(_pm_srcPath)+'/'+str(appName)

    dateTimeNow = time.strftime("%d-%m-%Y_%H.%M.%S")
    tuningReportFileName = str(_pm_srcPath)+'tuningReport'+'_'+str(dateTimeNow)+'.txt'
    tuningReportFile = open(tuningReportFileName, 'w')

    printParams(tuningReportFile)

    _pm_maxTime = 1000000


    # TODO: iterate over thread count for tuning                ( )

    # set the thread-count
    os.environ["OMP_NUM_THREADS"] = str(_pm_omp_threads)

    tuningTime_t1 = clock()

    globalMinConfig = 0
    globalMinTime = _pm_maxTime
    for _pm_config in range(1, _pm_configsCount+1):
        printLine()
        print "config #"+str(_pm_config),": ",

        # log to file
        printLine(tuningReportFile)
        print >>tuningReportFile, "Config #"+str(_pm_config),": ",

        # load shared library, name the function
        soFileName = str(progPrefix)+str(_pm_config)+'.so'
        libPipeline = ctypes.cdll.LoadLibrary(soFileName)
        pipelineFunc = libPipeline[libFunctionName]

        # TODO:
        # 1. Switch between report types - whether
        #    to print min / average / gmean etc.                ( )
        # 2. Graphics plots                                     ( )
        #

        # TODO:
        # 1. Use compiler info to get the correct prototype
        #    of the pipeline function to be called.             (X)
        # 2. Handle errors and dump them to log file            ( )

        # start timer
        localMinTime = _pm_maxTime
        for run in range(1, _pm_nruns+1):
            t1 = clock()
            pipelineFunc(*args)
            t2 = clock()
            t_total = float(t2) - float(t1)
            #print "time taken: ", t2*1000 - t1*1000, "ms"

            # local minima
            if float(t_total) < float(localMinTime):
                localMinTime = t_total
 
        # global minima
        if float(localMinTime) < float(globalMinTime):
            globalMinTime = localMinTime
            globalMinConfig = _pm_config

        print localMinTime*1000, \
              "(", globalMinTime*1000, ")"
        # write the same to log file
        # FIXME: fuse this with stdout dump, to work like 'tee'         ( )
        print >>tuningReportFile, localMinTime*1000

        # FIXME:
        # Unloading the shared lib is almost impossible. Nevertheless, try to
        # dispose off the loaded lib object
        del pipelineFunc
        _ctypes.dlclose(libPipeline._handle)
        del libPipeline

    tuningTime_t2 = clock()
    tuningTime = float(tuningTime_t2) - float(tuningTime_t1)

    printLine()
    print "Best Config :"
    print "Config #"+str(globalMinConfig)," -- ",
    print globalMinTime*1000
    print "Src Path    : \""+_pm_srcPath+"\""

    printLine()
    print "Tuning Time :"
    print tuningTime*1000

    printLine()

    # write the same to log file
    # FIXME: fuse this with stdout dump, to work like 'tee'             ( )
    printLine(tuningReportFile)
    print >>tuningReportFile, "Best Config :"
    print >>tuningReportFile, "Config #"+str(globalMinConfig)," -- ",
    print >>tuningReportFile, globalMinTime*1000
    print >>tuningReportFile, "Src Path    : \""+_pm_srcPath+"\""

    printLine(tuningReportFile)
    print >>tuningReportFile, "Tuning Time :"
    print >>tuningReportFile, tuningTime*1000

    printLine(tuningReportFile)
    tuningReportFile.close()

    return
