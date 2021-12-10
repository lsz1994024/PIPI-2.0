package proteomics.TheoSeq;

import proteomics.Types.SparseBooleanVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private final Map<String, Double> massTable = new HashMap<>();
    private final int missedCleavage;
    private double ms2Tolerance = 1.0005f;
    private double oneMinusBinOffset = 0.4f;
    private Pattern digestSitePattern;
    private final boolean cleavageFromCTerm;

    public MassTool(final int missedCleavage, Map<String, Double> fixModMap, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, double ms2Tolerance, double oneMinusBinOffset) {
        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm = cleavageFromCTerm;
        massTable.put("G", 57.021464f + fixModMap.get("G"));
        massTable.put("A", 71.037114f + fixModMap.get("A"));
        massTable.put("S", 87.032028f + fixModMap.get("S"));
        massTable.put("P", 97.052764f + fixModMap.get("P"));
        massTable.put("V", 99.068414f + fixModMap.get("V"));
        massTable.put("T", 101.047678f + fixModMap.get("I"));
        massTable.put("C", 103.009184f + fixModMap.get("C"));
        massTable.put("I", 113.084064f + fixModMap.get("I"));
        massTable.put("L", 113.084064f + fixModMap.get("L"));
        massTable.put("N", 114.042927f + fixModMap.get("N"));
        massTable.put("D", 115.026943f + fixModMap.get("D"));
        massTable.put("Q", 128.058578f + fixModMap.get("Q"));
        massTable.put("K", 128.094963f + fixModMap.get("K"));
        massTable.put("E", 129.042593f + fixModMap.get("E"));
        massTable.put("M", 131.040485f + fixModMap.get("M"));
        massTable.put("H", 137.058912f + fixModMap.get("H"));
        massTable.put("F", 147.068414f + fixModMap.get("F"));
        massTable.put("R", 156.101111f + fixModMap.get("R"));
        massTable.put("Y", 163.063329f + fixModMap.get("Y"));
        massTable.put("W", 186.079313f + fixModMap.get("W"));
        massTable.put("U", 150.953636f + fixModMap.get("U"));
        massTable.put("O", 237.147727f + fixModMap.get("O"));
        massTable.put("n", fixModMap.get("n"));
        massTable.put("#", 113.084064); // for I and L.
        massTable.put("$", 128.0767705); // for Q and K.
        massTable.put("C13_DIFF", 1.00335483);
        massTable.put("H2O", 18.010564684);
        massTable.put("NH3", 17.026549106);
        massTable.put("PROTON", 1.00727646688);
        massTable.put("Hatom", 1.007825032);
        massTable.put("Natom", 14.00307401);
        massTable.put("Oatom", 15.99491462);
        massTable.put("Patom", 30.97376151);
        massTable.put("Satom", 31.97207069);

        if (protectionSite.contentEquals("-")) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
        } else if (cleavageFromCTerm) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "](?=[^" + protectionSite + "])");
        } else {
            digestSitePattern = Pattern.compile("(?<=[^" + protectionSite + "])" + "[" + cleavageSite + "]");
        }
    }

    public double calResidueMass(String seq){ // Caution: nterm mod is not included!!!!!
        double totalMass = 0;
        int length = seq.length();
        for (int idx = 0; idx < length; ++idx) {
            totalMass += massTable.get(seq.substring(idx, idx + 1));
        }

        return totalMass;
    }

    public Set<String> buildPeptideSet(String proSeq) {
        Map<Integer, List<int[]>> digestRangeMap = digestTrypsin(proSeq);
        Set<String> peptideSeqSet = new HashSet<>();

        for (int i = 0; i <= missedCleavage; ++i) {
            for (int[] digestRange1 : digestRangeMap.get(i)) {
                String subString = proSeq.substring(digestRange1[0], digestRange1[1]);
                peptideSeqSet.add(subString);
            }
        }
        return peptideSeqSet;
    }

    public double[][] buildIonArray(String pepPeptide, int maxCharge) {
        // [NOTE] The b/y-ions charge 0
        double[][] peptideIonArray = new double[2 * maxCharge][pepPeptide.length()];
        double bIonMass = massTable.get("n");
        double yIonMass = calResidueMass(pepPeptide) + massTable.get("n") + massTable.get("H2O");

        for (int charge = 1; charge <= maxCharge; ++charge) {
            double bIonMassCharge = (bIonMass + charge * massTable.get("PROTON")) / charge;
            double yIonMassCharge = (yIonMass + charge * massTable.get("PROTON")) / charge;
            int chargeIdx = charge - 1;
            int chargeIdx2 = 2 * chargeIdx;
            int chargeIdx21 = chargeIdx2 + 1;

            for (int idx = 0; idx < pepPeptide.length(); ++idx) {
                // y-ion
                peptideIonArray[chargeIdx21][idx] = yIonMassCharge;

                String aa = pepPeptide.substring(idx, idx + 1);

                // b-ion
                bIonMassCharge += massTable.get(aa) / charge;
                peptideIonArray[chargeIdx2][idx] = bIonMassCharge;

                // Calculate next y-ion:
                if (idx == 0) {
                    yIonMassCharge -= (massTable.get(aa) / charge + massTable.get("n") / charge);
                } else {
                    yIonMassCharge -= massTable.get(aa) / charge;
                }
            }
        }

        return peptideIonArray;
    }

    public SparseBooleanVector buildVector(double[][] ionMatrix, int precursorCharge) {
        int colNum = ionMatrix[0].length;
        int rowNum = Math.min(ionMatrix.length / 2, precursorCharge - 1) * 2;
        if (precursorCharge == 1) {
            rowNum = 2;
        }

        List<Integer> idxList = new LinkedList<>();
        int maxIdx = 0;
        for (int i = 0; i < rowNum; ++i) {
            for (int j = 0; j < colNum; ++j) {
                if (ionMatrix[i][j] > 1e-6) {
                    int idx = mzToBin(ionMatrix[i][j]);
                    idxList.add(idx);
                    if (idx > maxIdx) {
                        maxIdx = idx;
                    }
                }
            }
        }

        SparseBooleanVector theoIonVector = new SparseBooleanVector();
        for (int idx : idxList) {
            theoIonVector.put(idx);
        }

        return theoIonVector;
    }

    public Map<String, Double> returnMassTable() {
        return massTable;
    }

    public int mzToBin(double mz) {
        return (int) Math.floor(mz / ms2Tolerance + oneMinusBinOffset);
    }

    public double binToMz(int idx) {
        return (idx - oneMinusBinOffset) * ms2Tolerance;
    }

    private Map<Integer, List<int[]>> digestTrypsin(String proSeq) {
        // Cut a protein
        List<Integer> cutPointList = new ArrayList<>(200);
        int length = proSeq.length();
        int idxStart = 0;
        Matcher matchObj = digestSitePattern.matcher(proSeq);
        cutPointList.add(0);
        while (idxStart < length) {
            if (matchObj.find()) {
                int cutPoint;
                if (cleavageFromCTerm) {
                    cutPoint = matchObj.end();
                } else {
                    cutPoint = matchObj.start();
                }
                cutPointList.add(cutPoint);
                idxStart = cutPoint;
            } else {
                cutPointList.add(length);
                break;
            }
        }

        Collections.sort(cutPointList);

        // Deal with missed cleavage
        Map<Integer, List<int[]>> digestRangeMap = new HashMap<>();
        for (int time = 0; time <= missedCleavage; ++time) {
            List<int[]> temp = new LinkedList<>();
            int leftPoint;
            int rightPoint;
            for (int i = 0; i + 1 + time < cutPointList.size(); ++i) {
                leftPoint = cutPointList.get(i);
                rightPoint = cutPointList.get(i + 1 + time);
                temp.add(new int[]{leftPoint, rightPoint});
            }
            digestRangeMap.put(time, temp);
        }

        return digestRangeMap;
    }
}
