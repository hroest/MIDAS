#!/usr/bin/python 
import sys, getopt, warnings, os, re

def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:o:s:p:e:d:",
                                    ["help",
                                     "working-dir",
                                     "output-file",
                                     "score_threshold",
                                     "spectral_count_threshold",
                                     "explain_peaks_threshold",
                                     "delta_threshold"])


    # Default working dir and config file
    working_dir = "./"
    outputFileName = ""
    dScore_Threshold = -1
    iSpectralCount_Threshold = -1
    iExplainedPeak_Threshold = -1
    dDeltaZ_Threshold = 0.02
    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory -o output_file_name (optional) -s score_threshold -p spectral_count_threshold -e explain_peaks_threshold -d delta_threshold"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-file"):
            outputFileName = value
        if option in ("-s", "--score_threshold"):
            dScore_Threshold = float(value) 
        if option in ("-p", "--spectral_count_threshold"):
            iSpectralCount_Threshold = int(value)
        if option in ("-e", "--explain_peaks_threshold"):
            iExplainedPeak_Threshold = int(value)
        if option in ("-d", "--delta_threshold"):
            #print value
            dDeltaZ_Threshold = float(value)

    # only -w is provided
    if (outputFileName == "") :
        outputFileName = working_dir+"goodpicks.txt"
    if ((dScore_Threshold <0) or (iSpectralCount_Threshold<0) or (iExplainedPeak_Threshold<0) or (dDeltaZ_Threshold < 0)) :
        print "please specify thresholds"
        sys.exit(1)
    AFT2_filename_list = get_file_list_with_ext(working_dir, ".AFT2")

    return outputFileName, AFT2_filename_list, dScore_Threshold, iSpectralCount_Threshold, iExplainedPeak_Threshold, dDeltaZ_Threshold



## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):
            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list


def getColumnId(sColumnNameLine, ColumnName) :
# This function returns column id (starting with 0) according to column name
# sColumnNameLine: the comlumn title line; ColumnName: column name
	ColumnName_list = sColumnNameLine.split("\t")
	try:
		iColumnId = ColumnName_list.index(ColumnName)
	except ValueError:
		print "can't find column "+ColumnName
		sys.exit(0)
	#print iColumnId
	return iColumnId

def ParseMebLine(sCurrentMebLine, iExplainedPeaks_ColumnId, iScore_ColumnId, iRank_ColumnId, iIdentifier_ColumnId, iFilename_ColumnId, iScanNumber_ColumnId) :
    Meb_info_list = sCurrentMebLine.split("\t")
    sExplainPeakInfo = Meb_info_list[iExplainedPeaks_ColumnId]
    sPeak_info_list = sExplainPeakInfo.split(">")
    iCurrentExplainedPeaks = int(sPeak_info_list[0])
    iCurrentTotalPeaks = int(sPeak_info_list[1])
    dCurrentScore = float(Meb_info_list[iScore_ColumnId])
    sCurrentIdentifier = Meb_info_list[iIdentifier_ColumnId]
    iCurrentRank  = int(Meb_info_list[iRank_ColumnId])
    sCurrentFilename = Meb_info_list[iFilename_ColumnId]
    sCurrentScanNumber = Meb_info_list[iScanNumber_ColumnId]
    return dCurrentScore, sCurrentIdentifier, iCurrentExplainedPeaks, iCurrentTotalPeaks, iCurrentRank, sCurrentFilename, sCurrentScanNumber

def ParseAFT2Line(sCurrentAFT2Line, iMofZ_ColumnId, iIntensity_ColumnId, iSmiles_ColumnId) :
    match_info_list = sCurrentAFT2Line.split("\t")
    dCurrentMofZ = float(match_info_list[iMofZ_ColumnId])
    dCurrentIntensity = float(match_info_list[iIntensity_ColumnId])
    sCurrentSmiles = match_info_list[iSmiles_ColumnId]
    return dCurrentMofZ, dCurrentIntensity, sCurrentSmiles 

def ReadAFT2Files(AFT2_filename_list, dScore_Threshold, iExplainedPeak_Threshold, dDeltaZ_Threshold):
    dScoreThreshold = dScore_Threshold#0.2
    iExplainedPeakThreshold  = iExplainedPeak_Threshold  #5
    dDeltaZThreshold= dDeltaZ_Threshold  #0.02
    Pick_dict = {}
    for each_filename in AFT2_filename_list:
        iLineNum = 0
        bFirstEntry = True
        current_AFT2_file = open(each_filename)
        sCurrent_AFT2_list = []
        for each_line in current_AFT2_file :
            each_line = each_line.strip()
            if (each_line == "") :
                continue
            iLineNum += 1
            if (iLineNum == 1) :
                iExplainedPeaks_ColumnId = getColumnId(each_line, "ExplainedPeaks") + 1
                iScore_ColumnId = getColumnId(each_line, "Score") + 1
                iRank_ColumnId  = getColumnId(each_line, "Rank") + 1
                iIdentifier_ColumnId  = getColumnId(each_line, "Identifier") + 1
                iFilename_ColumnId    = getColumnId(each_line, "Filename") + 1
                iScanNumber_ColumnId  = getColumnId(each_line, "ScanNumber") + 1
                sFirstTitleLine = each_line
            elif (iLineNum == 2) :
                iMofZ_ColumnId = getColumnId(each_line, "m/z")
                iIntensity_ColumnId = getColumnId(each_line, "Intensity")
                iSmiles_ColumnId = getColumnId(each_line, "SMILES")
                sSecondTitleLine = each_line
            elif (each_line.startswith("M\t")) :
                if not (bFirstEntry) :
                    sPreviousMebLine    = sCurrentMebLine
                    dPreviousScore      = dCurrentScore
                    sPreviousIdentifier = sCurrentIdentifier
                    iPreviousExplainedPeaks = iCurrentExplainedPeaks
                    iPreviousTotalPeaks = iCurrentTotalPeaks
                    sPrevious_AFT2_list = sCurrent_AFT2_list
                    iPreviousRank       = iCurrentRank
                    dPreviousDeltaZ     = dPreviousScore
                    sPreviousFilename   = sCurrentFilename
                    sPreviousScanNumber = sCurrentScanNumber
                sCurrentMebLine = each_line
                dCurrentScore, sCurrentIdentifier, iCurrentExplainedPeaks, iCurrentTotalPeaks, iCurrentRank, sCurrentFilename, sCurrentScanNumber = ParseMebLine(sCurrentMebLine, iExplainedPeaks_ColumnId, iScore_ColumnId, iRank_ColumnId, iIdentifier_ColumnId, iFilename_ColumnId, iScanNumber_ColumnId)
                if not (bFirstEntry) :
                    if (sPreviousFilename == sCurrentFilename) and (sPreviousScanNumber == sCurrentScanNumber) :
                        dPreviousDeltaZ = dPreviousScore - dCurrentScore
                        #print sPreviousRank 
                        #if (sPreviousRank == "1") :
                        #    print dPreviousScore, dPreviousDeltaZ
                    #print iPreviousRank, dPreviousScore, iPreviousExplainedPeaks, dPreviousDeltaZ
                    if ((iPreviousRank == 1) and (dPreviousScore >= dScoreThreshold) and (iPreviousExplainedPeaks >= iExplainedPeakThreshold) and (dPreviousDeltaZ >= dDeltaZThreshold)) :
            #            print dDeltaZ_Threshold
                        if sPreviousIdentifier in Pick_dict :
                            exist_entry = Pick_dict[sPreviousIdentifier]
                            dExistScore = exist_entry[2]
                            dExistSpectralCount = exist_entry[3]
                            if (dPreviousScore>dExistScore) :
                                Pick_dict[sPreviousIdentifier] = [sPreviousMebLine, sPrevious_AFT2_list, dPreviousScore, dExistSpectralCount+1]
                            else :
                                exist_entry[3] += 1
                        else :
                            Pick_dict[sPreviousIdentifier] = [sPreviousMebLine, sPrevious_AFT2_list, dPreviousScore, 1]#3rd element is spectral count
                else :
                    bFirstEntry = False
                sCurrent_AFT2_list = []
                sCurrentRank = "-1"
            else:
                sCurrentAFT2Line = each_line
                dCurrentMofZ, dCurrentIntensity, sCurrentSmiles = ParseAFT2Line(sCurrentAFT2Line, iMofZ_ColumnId, iIntensity_ColumnId, iSmiles_ColumnId)
                sCurrent_AFT2_list.append([sCurrentAFT2Line, dCurrentMofZ, dCurrentIntensity, sCurrentSmiles])
        if not (bFirstEntry) :
            sPreviousMebLine    = sCurrentMebLine
            dPreviousScore      = dCurrentScore
            sPreviousIdentifier = sCurrentIdentifier
            iPreviousExplainedPeaks = iCurrentExplainedPeaks
            iPreviousTotalPeaks = iCurrentTotalPeaks
            sPrevious_AFT2_list = sCurrent_AFT2_list
            iPreviousRank       = iCurrentRank
            dPreviousDeltaZ     = dPreviousScore
            sPreviousFilename   = sCurrentFilename
            sPreviousScanNumber = sCurrentScanNumber
            if ((iPreviousRank == 1) and (dPreviousScore >= dScoreThreshold) and (iPreviousExplainedPeaks >= iExplainedPeakThreshold) and (dPreviousDeltaZ >= dDeltaZThreshold)) :
                if sPreviousIdentifier in Pick_dict :
                    exist_entry = Pick_dict[sPreviousIdentifier] 
                    dExistScore = exist_entry[2]
                    dExistSpectralCount = exist_entry[3]
                    if (dPreviousScore>dExistScore) :
                        Pick_dict[sPreviousIdentifier] = [sPreviousMebLine, sPrevious_AFT2_list, dPreviousScore, dExistSpectralCount+1]
                    else :
                        exist_entry[3] += 1
                else :
                    Pick_dict[sPreviousIdentifier] = [sPreviousMebLine, sPrevious_AFT2_list, dPreviousScore, 1]
        current_AFT2_file.close()
        return Pick_dict, sFirstTitleLine, sSecondTitleLine

def Filtering(Pick_dict, outputFileName, sFirstTitleLine, sSecondTitleLine, iSpectralCount_Threshold) :
    iSpectralCountThreshold = iSpectralCount_Threshold  #3
    output_file = open(outputFileName, "w")
    output_file.write(sFirstTitleLine+"\n")
    output_file.write(sSecondTitleLine+"\n")
    for sCurrentIdentifier, CurrentValue_list in Pick_dict.iteritems() :
        iCurrentSpectralCount = CurrentValue_list[3]
        #print iCurrentSpectralCount
        if (iCurrentSpectralCount >= iSpectralCountThreshold) :
            sCurrentMebLine = CurrentValue_list[0]
            output_file.write(sCurrentMebLine+"\n")
            currentAFT2_list = CurrentValue_list[1]
            for each_peak_list in currentAFT2_list :
                sCurrentAFT2Line = each_peak_list[0]
                dCurrentMofZ     = each_peak_list[1]
                dCurrentIntensity= each_peak_list[2]
                sCurrentSmiles   = each_peak_list[3]
                #output_file.write(str(dCurrentMofZ)+"\t"+str(dCurrentIntensity)+"\t"+sCurrentSmiles+"\n")
                output_file.write(sCurrentAFT2Line+"\n")
    output_file.close()

def main(argv=None):
    if argv is None:
        argv = sys.argv
        outputFileName, AFT2_filename_list, dScore_Threshold, iSpectralCount_Threshold, iExplainedPeak_Threshold, dDeltaZ_Threshold = parse_options(argv)
        Pick_dict, sFirstTitleLine, sSecondTitleLine = ReadAFT2Files(AFT2_filename_list, dScore_Threshold, iExplainedPeak_Threshold, dDeltaZ_Threshold)
        #print len(Pick_dict)
        Filtering(Pick_dict, outputFileName, sFirstTitleLine, sSecondTitleLine, iSpectralCount_Threshold)

if __name__ == "__main__":
    main()
