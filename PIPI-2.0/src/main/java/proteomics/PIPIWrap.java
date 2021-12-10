package proteomics;

import proteomics.Index.BuildIndex;
import proteomics.Search.CalEValue;
import proteomics.Library.PrepareSpectrum;
import proteomics.Search.CalSubscores;
//import proteomics.Search.CalXcorr;
import proteomics.PTM.FindPTM;
import proteomics.Search.Search;
import proteomics.Library.Binomial;
import proteomics.Segment.InferenceSegment;
import proteomics.Spectrum.PreSpectra;
import proteomics.Library.MassTool;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<PIPIWrap.Entry> {
    private final InferenceSegment inference3SegmentObj;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double ms2Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    //    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double precursorMass;
    private final PrepareSpectrum preSpectrum;
    private final String sqlPath;
    private final Binomial binomial;


    public PIPIWrap(InferenceSegment inference3SegmentObj, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass, FindPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, Binomial binomial) {
        this.inference3SegmentObj = inference3SegmentObj;
        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.preSpectrum = preSpectrum;
        this.sqlPath = sqlPath;
        this.binomial = binomial;
//        peptide0Map = buildIndex.getPeptide0Map();
    }


    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }

        // preprocess peak list
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);  // remove reporter ions and keep the local top peaks

        if (plMap.isEmpty()) {
            return null;
        }
        // Coding
        InferenceSegment inferSegment = buildIndex.getInferSegment();
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(precursorMass, plMap);
        List<TwoExpAA> expAaLists2 = inference3SegmentObj.inferSegmentLocationFromSpectrum2(precursorMass, plMap);
        List<FourExpAA> expAaLists4 = inference3SegmentObj.inferSegmentLocationFromSpectrum4(precursorMass, plMap);


        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inferSegment.generateSegmentIntensityVector(expAaLists);
            SparseVector scanCode2 = inferSegment.generateSegmentIntensityVector2(expAaLists2);
            SparseVector scanCode4 = inferSegment.generateSegmentIntensityVector4(expAaLists4);

            // Begin search.
//            Search searchObj = new Search(buildIndex, precursorMass, scanCode, scanCode2, scanCode4, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);
            Search search = new Search(buildIndex, precursorMass, scanCode, scanCode2, scanCode4, inference3SegmentObj, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);
// prepare the spectrum
            SparseBooleanVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
            } else {
                expProcessedPL = preSpectrum.digitizePL(plMap);
            }

            // dynamic programming
            double localMS1ToleranceL = -1 * ms1Tolerance;
            double localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass * leftInverseMs1Tolerance) - precursorMass;
                localMS1ToleranceR = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
            }

            // infer PTM using the new approach


            if (!search.getPTMOnlyResult().isEmpty()) {
//                Peptide[] peptideArray = peptideSet.toArray(new Peptide[0]);
//                Peptide topPeptide = peptideArray[0];
//                TreeSet<Peptide> ptmPatterns = null;
//                if (topPeptide.hasVarPTM()) {
//                    ptmPatterns = modSequences.get(topPeptide.getPTMFreePeptide());
//                }
//                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge, ptmPatterns, binomial);

                Connection sqlConnection = DriverManager.getConnection(sqlPath);
                Statement sqlStatement = sqlConnection.createStatement();
                ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, isDecoy, score FROM spectraTable WHERE scanId='%s'", scanId));
                if (sqlResultSet.next()) {
                    int scanNum = sqlResultSet.getInt("scanNum");
                    String scanId = sqlResultSet.getString("scanId");
                    int precursorCharge = sqlResultSet.getInt("precursorCharge");
                    double precursorMass = sqlResultSet.getDouble("precursorMass");
                    String mgfTitle = sqlResultSet.getString("mgfTitle");
                    int isotopeCorrectionNum = sqlResultSet.getInt("isotopeCorrectionNum");
                    double ms1PearsonCorrelationCoefficient = sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient");

                    double deltaLCn = 1;
//                    if (peptideArray.length > 4) {
//                        deltaLCn = (peptideArray[0].getScore() - peptideArray[4].getScore()) / peptideArray[0].getScore();
//                    }
//                    double deltaCn = 1;
//                    if (peptideArray.length > 1) {
//                        deltaCn = (peptideArray[0].getScore() - peptideArray[1].getScore()) / peptideArray[0].getScore();
                }

                String otherPtmPatterns = "-";
//                    if (ptmPatterns != null) {
//                        List<String> tempList = new LinkedList<>();
//                        Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
//                        ptmPatternsIterator.next();
//                        while (ptmPatternsIterator.hasNext()) {
//                            Peptide temp = ptmPatternsIterator.next();
//                            tempList.add(String.format(Locale.US, "%s-%.4f", temp.getPtmContainingSeq(buildIndex.returnFixModMap()), temp.getScore())); // Using 4 decimal here because it is write the the result file for checking. It is not used in scoring or other purpose.
//                        }
//                        otherPtmPatterns = String.join(";", tempList);
//                    }

                String PTMFreeCandidates = "";
                if (search.getPTMFreeResult() != null) {
                    List<String> tempList = new LinkedList<>();
                    for (Peptide temp : search.getPTMFreeResult()) {
                        tempList.add(String.format(Locale.US, "%s", temp.getPTMFreeSeq()));
                    }
                    PTMFreeCandidates = String.join(";", tempList);
                }

                String PTMOnlyCandidates = "";
                if (search.getPTMOnlyResult() != null) {
                    List<String> tempList = new LinkedList<>();
                    for (Peptide temp : search.getPTMOnlyResult()) {
                        tempList.add(String.format(Locale.US, "%s", temp.getPTMFreeSeq()));
                    }
                    PTMOnlyCandidates = String.join(";", tempList);
                }

                Entry entry = new Entry(0, scanId, precursorCharge, precursorMass, "", 0, 0, "", search.getPTMOnlyResult().get(0).getPTMContainedString(buildIndex.returnFixModMap()), 0, search.getPTMOnlyResult().get(0).isDecoy() ? 1 : 0, 0,0, 0, 0.0, 0.0, 0, 0, 0, 0, otherPtmPatterns, "");

//                Entry entry = new Entry(scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, topPeptide.getPTMContainedString(buildIndex.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLCn, deltaCn, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), otherPtmPatterns, topPeptide.getaScore(), PTMFreeCandidates, PTMOnlyCandidates);

                sqlResultSet.close();
                sqlStatement.close();
                return entry;
            } else {
                throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
            }
        }
        return null;
    }


    public class Entry {

        final int scanNum;
        final String scanId;
        final int precursorCharge;
        final double precursorMass;
        final String mgfTitle;
        final int isotopeCorrectionNum;
        final double ms1PearsonCorrelationCoefficient;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        final int globalRank;
        final double normalizedCorrelationCoefficient;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns; // It has 4 decimal because it is write the the result file for checking. It is not used in scoring or other purpose.
        final String aScore;

        Entry(int scanNum, String scanId, int precursorCharge, double precursorMass, String mgfTitle, int isotopeCorrectionNum, double ms1PearsonCorrelationCoefficient, String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient, double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac, double explainedAaFrac, String otherPtmPatterns, String aScore) {
            this.scanNum = scanNum;
            this.scanId = scanId;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.mgfTitle = mgfTitle;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.globalRank = globalRank;
            this.normalizedCorrelationCoefficient = normalizedCorrelationCoefficient;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
        }
    }
}


