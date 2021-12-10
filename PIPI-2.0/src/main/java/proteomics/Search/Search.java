package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.Index.BuildIndex;
import proteomics.Library.MassTool;
import proteomics.PIPI;
import proteomics.Types.Peptide;
import proteomics.Types.ResultEntry;
import proteomics.Types.SparseBooleanVector;
import proteomics.Types.SparseVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Search {

    private static final Logger logger = LoggerFactory.getLogger(Search.class);
    private static final int rankNum = 10;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();
    private final MassTool massToolObj;
    private final int maxMs2Charge;

    public Search(BuildIndex buildIndexObj, double precursorMass, SparseVector numCode, SparseVector numCode2, SparseVector numCode4, InferenceSegment inference3SegmentObj, double ms1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int maxMs2Charge) {
        massToolObj = buildIndexObj.returnMassToolObj();
        this.maxMs2Charge = maxMs2Charge;

        Map<Character, Double> massTable = massToolObj.getMassTable();
        Set<String> targetPeptideSet = buildIndexObj.returnPepMassMap().keySet();
        Map<String, Double> localPeptideMassMap = new HashMap<>();
        localPeptideMassMap.putAll(buildIndexObj.returnPepMassMap());
        localPeptideMassMap.putAll(buildIndexObj.returnDecoyPepMassMap());
        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> proteinSeqMap = buildIndexObj.returnProPepMap();

        Map<Integer, Double> numCodeNormSquareMap = new HashMap<>();
        Map<Integer, Double> numCodeNormSquareMap2 = new HashMap<>();
        Map<Integer, Double> numCodeNormSquareMap4 = new HashMap<>();

        SparseVector code = numCode;
        SparseVector code2 = numCode2;
        SparseVector code4 = numCode4;
        double threeNorm = code.norm2square();
        double twoNorm = code2.norm2square();
        double fourNorm = code4.norm2square();

        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);


        for (String peptide : localPeptideMassMap.keySet()) {
            double peptideMass = localPeptideMassMap.get(peptide);
            SparseBooleanVector peptideThreeCode = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(peptide));
            SparseBooleanVector peptideTwoCode = inference3SegmentObj.generateSegmentBooleanVector2(inference3SegmentObj.cutTheoSegment2(peptide));
            SparseBooleanVector peptideFourCode = inference3SegmentObj.generateSegmentBooleanVector4(inference3SegmentObj.cutTheoSegment4(peptide));
            double peptideCodeThreeNormSquare = peptideThreeCode.norm2square();
            double peptideCodeTwoNormSquare = peptideTwoCode.norm2square();
            double peptideCodeFourNormSquare = peptideFourCode.norm2square();

            double tolerance = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                tolerance = peptideMass * ms1Tolerance / 1e6f;
            }

            double leftMass = peptideMass + minPtmMass;
            double rightMass = peptideMass + maxPtmMass;

            // The peptide mass is too small or too large.
            if (leftMass >= rightMass) {
                continue;
            }

            boolean isTarget = false;
            if (targetPeptideSet.contains(peptide)) {
                isTarget = true;
            }
            double score = 0.0;
            double score2 = 0.0;
            double score4 = 0.0;
            double overallscore = 0.0;

            double temp1 = Math.sqrt(peptideCodeThreeNormSquare * threeNorm);
            double temp2 = Math.sqrt(peptideCodeTwoNormSquare * twoNorm);
            double temp4 = Math.sqrt(peptideCodeFourNormSquare * fourNorm);
//                        if (temp1 > 1e-6) {
//                            score = peptideThreeCode.dot(scanCode) / temp1;
////                            score2 = peptideTwoCode.dot(scanCode2) / temp2;
//                            score4 = peptideFourCode.dot(scanCode4) / temp4;
//                        }
            if (temp1 > 1e-6) {
                score = peptideThreeCode.dot(code) / temp1;
            }
            if (temp2 > 1e-6) {
                score2 = peptideTwoCode.dot(code2) / temp2;
            }
            if (temp4 > 1e-6) {
                score4 = peptideFourCode.dot(code4) / temp4;
            }

//                        overallscore = score4;
            overallscore = score + score2 + score4;
            double deltaMass = precursorMass - peptideMass;


            if (isTarget) {
                if (Math.abs(deltaMass) <= tolerance) {
                                // PTM-free
                    if (ptmFreeQueue.size() < rankNum) {
                        ptmFreeQueue.add(new ResultEntry(overallscore, peptide, false, false));
                    }
                    else {
                        if (overallscore > ptmFreeQueue.peek().score) {
                            ptmFreeQueue.poll();
                            ptmFreeQueue.add(new ResultEntry(overallscore, peptide, false, false));
                        }
                    }
                }

                if (Math.abs(deltaMass) > tolerance) {
                                // PTM-only
                    if (ptmOnlyQueue.size() < rankNum) {
                        ptmOnlyQueue.add(new ResultEntry(overallscore, peptide, false, true));
                    } else {
                        if (overallscore > ptmOnlyQueue.peek().score) {
                            ptmOnlyQueue.poll();
                            ptmOnlyQueue.add(new ResultEntry(overallscore, peptide, false, true));
                        }
                    }
                }
            } else {
                if (Math.abs(deltaMass) <= tolerance) {
                                // PTM-free
                    if (ptmFreeQueue.size() < rankNum) {
                        ptmFreeQueue.add(new ResultEntry(overallscore, peptide, true, false));
                    } else {
                        if (overallscore > ptmFreeQueue.peek().score) {
                            ptmFreeQueue.poll();
                            ptmFreeQueue.add(new ResultEntry(overallscore, peptide, true, false));
                        }
                    }
                }

                if (Math.abs(deltaMass) > tolerance) {
                                // PTM-only
                    if (ptmOnlyQueue.size() < rankNum) {
                        ptmOnlyQueue.add(new ResultEntry(overallscore, peptide, true, true));
                    } else {
                        if (overallscore > ptmOnlyQueue.peek().score) {
                            ptmOnlyQueue.poll();
                            ptmOnlyQueue.add(new ResultEntry(overallscore, peptide, true, true));
                        }
                    }
                }
            }
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massToolObj, maxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massToolObj, maxMs2Charge);
        }
    }

//    private void writeResult(String filePath,  Map<Integer, List<Peptide>> ptmOnly) {
//        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
//            writer.write("scanNum,peptide,globalRank,is_decoy\n");
//            for (int scanNum : ptmOnly.keySet()) {
//                List<Peptide> peptideList = ptmOnly.get(scanNum);
//                for (Peptide peptide : peptideList) {
//                    if (peptide.isDecoy()) {
//                        writer.write(scanNum + "," + peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",1\n");
//                    } else {
//                        writer.write(scanNum + "," + peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",0\n");
//                    }
//                }
//            }
//        } catch (IOException ex) {
//            ex.printStackTrace();
//            logger.error(ex.getMessage());
//            System.exit(1);
//        }
//    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, globalRank));
            --globalRank;
        }

        return peptideList;
    }

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }
}
