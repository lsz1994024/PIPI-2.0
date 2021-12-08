package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.FindPTM;
import proteomics.Search.CalEValue;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalXcorr;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;
import java.util.concurrent.Callable;

public class PIPIWrap implements Callable<List<FinalResultEntry>> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final InferenceSegment inference3SegmentObj;
//    private final InferenceSegment inference2SegmentObj;
//    private final InferenceSegment inference4SegmentObj;
    private final Map<Integer, SpectrumEntry> numSpectrumMap;
    private final NavigableMap<Float, List<Integer>> subMassNumMap;
    private final Map<String, TreeSet<Integer>> siteMass1000Map;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final int batchStartIdx;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, InferenceSegment inference3SegmentObj, Map<Integer, SpectrumEntry> numSpectrumMap, NavigableMap<Float, List<Integer>> subMassNumMap, Map<String, TreeSet<Integer>> siteMass1000Map, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, int batchStartIdx) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.inference3SegmentObj = inference3SegmentObj;
//        this.inference2SegmentObj = inference2SegmentObj;
//        this.inference4SegmentObj = inference4SegmentObj;
        this.numSpectrumMap = numSpectrumMap;
        this.subMassNumMap = subMassNumMap;
        this.siteMass1000Map = siteMass1000Map;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        this.batchStartIdx = batchStartIdx;
    }

    @Override
    public List<FinalResultEntry> call() {
        try {
            // Inference amino acids.
            Map<Integer, SparseVector> numCodeMap = new HashMap<>();
            Map<Integer, SparseVector> numCodeMap2 = new HashMap<>();
            Map<Integer, SparseVector> numCodeMap4 = new HashMap<>();
            Map<Integer, List<ThreeExpAA>> numExp3aaLists = new HashMap<>();
            Map<Integer, List<TwoExpAA>> numExp2aaLists = new HashMap<>();
            Map<Integer, List<FourExpAA>> numExp4aaLists = new HashMap<>();
            for (List<Integer> numList : subMassNumMap.values()) {
                for (int scanNum : numList) {
                    SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);

                    // Coding
                    List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(spectrumEntry);
                    List<TwoExpAA> expAaLists2 = inference3SegmentObj.inferSegmentLocationFromSpectrum2(spectrumEntry);
                    List<FourExpAA> expAaLists4 = inference3SegmentObj.inferSegmentLocationFromSpectrum4(spectrumEntry);


                    if (!expAaLists.isEmpty()) {
                        numExp3aaLists.put(scanNum, expAaLists);
                        numExp2aaLists.put(scanNum, expAaLists2);
                        numExp4aaLists.put(scanNum, expAaLists4);
                        numCodeMap.put(scanNum, inference3SegmentObj.generateSegmentIntensityVector(expAaLists));
                        numCodeMap2.put(scanNum, inference3SegmentObj.generateSegmentIntensityVector2(expAaLists2));
                        numCodeMap4.put(scanNum, inference3SegmentObj.generateSegmentIntensityVector4(expAaLists4));
                    }
                }
            }

            // Begin search.
            Search searchObj = new Search(buildIndexObj, numCodeMap, numCodeMap2, numCodeMap4, inference3SegmentObj, subMassNumMap, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge, batchStartIdx);

            Set<Integer> ScanNumSet = new HashSet<Integer>();
            Set<Integer> ptmFreeScanNumSet = searchObj.getPTMFreeResult().keySet();
            Set<Integer> ptmOnlyScanNumSet = searchObj.getPTMOnlyResult().keySet();
            ScanNumSet.addAll(ptmFreeScanNumSet);
            ScanNumSet.addAll(ptmOnlyScanNumSet);
            for (int ScanNum: ScanNumSet) {
                String ptmFreeSequence = null;
                double ptmFreeScore = 0.0;
                String ptmOnlySequence = null;
                double ptmOnlyScore = 0.0;
                if (searchObj.getPTMFreeResult().containsKey(ScanNum)) {
                    LinkedList<Peptide> ptmFreeList = (LinkedList<Peptide>) searchObj.getPTMFreeResult().get(ScanNum);
                    ptmFreeSequence = ptmFreeList.get(ptmFreeList.size()-1).getPTMFreeSeq();
                    ptmFreeScore = ptmFreeList.get(ptmFreeList.size()-1).getNormalizedCrossCorr();
                }
                if (searchObj.getPTMOnlyResult().containsKey(ScanNum)) {
                    LinkedList<Peptide> ptmOnlyList = (LinkedList<Peptide>) searchObj.getPTMOnlyResult().get(ScanNum);
                    ptmOnlySequence = ptmOnlyList.get(ptmOnlyList.size()-1).getPTMFreeSeq();
                    ptmOnlyScore = ptmOnlyList.get(ptmOnlyList.size()-1).getNormalizedCrossCorr();
                }
//                if (ptmFreeScore>=ptmOnlyScore) {
//                    System.out.println(ScanNum+" "+ptmFreeSequence);
//                }
                if (ptmFreeScore<ptmOnlyScore) {
//                    System.out.println(ScanNum+" "+ptmOnlySequence);
                    System.out.println(ScanNum+" "+ptmFreeSequence);
                }
            }

            logger.debug("Analyzing PTMs...");
//            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), numSpectrumMap, numExp3aaLists, massToolObj, siteMass1000Map, minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, batchStartIdx);
//            Map<Integer, List<Peptide>> ptmOnlyTemp = findPtmObj.getPeptidesWithPTMs();

            Map<Integer, List<Peptide>> ptmOnlyTemp = searchObj.getPTMOnlyResult();

            logger.debug("Calculating final score...");
            // put ptmFree and ptmOnly together
            Map<Integer, List<Peptide>> numCandidateMap = new HashMap<>(searchObj.getPTMFreeResult());
            for (int scanNum : ptmOnlyTemp.keySet()) {
                if (numCandidateMap.containsKey(scanNum)) {
                    numCandidateMap.get(scanNum).addAll(ptmOnlyTemp.get(scanNum));
                } else {
                    numCandidateMap.put(scanNum, ptmOnlyTemp.get(scanNum));
                }
            }
            CalXcorr calXcorrObj = new CalXcorr(numCandidateMap, numSpectrumMap, massToolObj, buildIndexObj);
            List<FinalResultEntry> subScoredPsms = calXcorrObj.getScoredPSMs();

            new CalSubscores(subScoredPsms, numSpectrumMap, ms2Tolerance);

            logger.debug("Estimating E-value for each PSM...");
            for (FinalResultEntry psm : subScoredPsms) {
                new CalEValue(psm);
            }
            return subScoredPsms;
        } catch (Exception ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }
        return new LinkedList<>();
    }
}
