package proteomics.Index;

import java.util.*;

import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.*;
import proteomics.Library.MassTool;

public class BuildIndex {

    private static final int maxLoopTime = 10;

    private double minPrecursorMass = 0;
    private double maxPrecursorMass = 0;
    private final MassTool massTool;
    private Map<Character, Double> massTable;
    private Map<String, String> proPeptideMap = new HashMap<>();
    private Map<String, Double> peptideMassMap = new HashMap<>();
    private InferenceSegment inferSegment;
    private Map<String, Set<String>> peptideProMap = new HashMap<>();
    private Set<String> forCheckDuplicate = new HashSet<>();
    private Map<String, Double> decoyPeptideMassMap = new HashMap<>();
    private Map<String, String> decoyPeptideProMap = new HashMap<>();
    private Map<Character, Double> fixModMap = new HashMap<>();
    private TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();

    /////////////////////////////////public methods//////////////////////////////////////////////////////////////////
    public BuildIndex(Map<String, String> parameterMap, String labelling) {
        // initialize parameters
        minPrecursorMass = Double.valueOf(parameterMap.get("min_precursor_mass"));
        maxPrecursorMass = Double.valueOf(parameterMap.get("max_precursor_mass"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));

        // Read fix modification
        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));

        // read protein database
        DbTool dbToolObj = new DbTool(dbPath);
        proPeptideMap = dbToolObj.returnSeqMap();

        // define a new MassTool object
        massTool = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site_1").trim(), parameterMap.get("protection_site_1").trim(), parameterMap.get("is_from_C_term_1").trim().contentEquals("1"), parameterMap.getOrDefault("cleavage_site_2", null), parameterMap.getOrDefault("protection_site_2", null), parameterMap.containsKey("is_from_C_term_2") ? parameterMap.get("is_from_C_term_2").trim().contentEquals("1") : null, ms2Tolerance, oneMinusBinOffset, labelling);
        massTable = massTool.getMassTable();

        // build database
        buildPeptideMap();
        buildDecoyPepChainMap();
        buildMassPeptideMap();
    }

    /////////////////////////////////////public methods////////////////////////////////////////////////////////////////////
    public MassTool returnMassToolObj() {
        return massTool;
    }

    public Map<String, String> returnProPepMap() {
        return proPeptideMap;
    }

    public Map<String, Double> returnPepMassMap() {
        return peptideMassMap;
    }

    public Map<String, Set<String>> returnPepProMap() {
        return peptideProMap;
    }

    public InferenceSegment getInferSegment() {
        return inferSegment;
    }
    public Map<String, Double> returnDecoyPepMassMap() {
        return decoyPeptideMassMap;
    }

    public Map<String, String> returnDecoyPepProMap() {
        return decoyPeptideProMap;
    }

    public Map<Character, Double> returnFixModMap() {
        return fixModMap;
    }

    public TreeMap<Double, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    private void buildPeptideMap() {
        Set<String> proIdSet = proPeptideMap.keySet();
        for (String proId : proIdSet) {
            String proSeq = proPeptideMap.get(proId);
            Set<String> peptideSet = massTool.buildPeptideSet(proSeq);
            for (String peptide : peptideSet) {
                if (peptide.contains("B") || peptide.contains("J") || peptide.contains("X") || peptide.contains("Z")) {
                    continue;
                }

                double massTemp = massTool.calResidueMass(peptide) + massTable.get("n") + massTable.get("H2O"); // calMass just calculate the residue mass, so we should add a H2O
                if ((massTemp <= maxPrecursorMass) && (massTemp >= minPrecursorMass)) {
                    peptideMassMap.put(peptide, massTemp);

                    // Add the sequence to the check set for decoy duplicate check
                    String templateSeq = peptide.replace('L', 'I'); // "L" and "I" have the same mass.
                    forCheckDuplicate.add(templateSeq);

                    if (peptideProMap.containsKey(peptide)) {
                        Set<String> proList = peptideProMap.get(peptide);
                        proList.add(proId);
                        peptideProMap.put(peptide, proList);
                    } else {
                        Set<String> proList = new HashSet<>();
                        proList.add(proId);
                        peptideProMap.put(peptide, proList);
                    }
                }
            }
        }
    }

    private void buildDecoyPepChainMap() {
        Set<String> peptideSet = peptideProMap.keySet();
        for (String originalPeptide : peptideSet) {
            String decoyPeptide = shuffleSeq2(originalPeptide);
            if (decoyPeptide.isEmpty()) {
                continue;
            }

            double decoyMassTemp = peptideMassMap.get(originalPeptide);
            decoyPeptideMassMap.put(decoyPeptide, decoyMassTemp);
            String proId = peptideProMap.get(originalPeptide).iterator().next();
            String decoyProId = "DECOY_" + proId;
            decoyPeptideProMap.put(decoyPeptide, decoyProId);
        }
    }

    private void buildMassPeptideMap() {
        for (String peptide : peptideMassMap.keySet()) {
            double mass = peptideMassMap.get(peptide);
            if (massPeptideMap.containsKey(mass)) {
                massPeptideMap.get(mass).add(peptide);
            } else {
                Set<String> temp = new HashSet<>();
                temp.add(peptide);
                massPeptideMap.put(mass, temp);
            }
        }
    }

    private String reverseSeq(String seq) {
        String decoySeq;
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            StringBuilder sb = new StringBuilder(seq.substring(0, seq.length() - 1)).reverse();
            decoySeq = sb.toString() + seq.substring(seq.length() - 1);
        } else {
            StringBuilder sb = new StringBuilder(seq).reverse();
            decoySeq = sb.toString();
        }
        if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
            return "";
        } else {
            return decoySeq;
        }
    }

    private String shuffleSeq(String seq) {
        List<Character> seqList = new ArrayList<>(seq.length() + 2);
        Random randomObj = new Random();
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            for (int i = 0; i < seq.length() - 1; ++i) {
                seqList.add(seq.charAt(i));
            }
            StringBuilder sb;
            int loopTime = maxLoopTime;
            do {
                --loopTime;
                Collections.shuffle(seqList, randomObj);
                sb = new StringBuilder(seq.length() + 2);
                for (char aa : seqList) {
                    sb.append(aa);
                }
                sb.append(seq.charAt(seq.length() - 1));
            } while ((forCheckDuplicate.contains(sb.toString().replace('L', 'I'))) && (loopTime > 0));
            return sb.toString();
        } else {
            for (int i = 0; i < seq.length(); ++i) {
                seqList.add(seq.charAt(i));
            }
            StringBuilder sb;
            int loopTime = maxLoopTime;
            do {
                --loopTime;
                Collections.shuffle(seqList, randomObj);
                sb = new StringBuilder(seq.length() + 2);
                for (char aa : seqList) {
                    sb.append(aa);
                }
            } while ((forCheckDuplicate.contains(sb.toString().replace('L', 'I'))) && (loopTime > 0));
            return sb.toString();
        }
    }

    private String shuffleSeq2(String seq) {
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            char[] tempArray = seq.substring(0, seq.length() - 1).toCharArray();
            int idx = 0;
            while (idx < tempArray.length - 1) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            }
            String decoySeq = String.valueOf(tempArray) + seq.substring(seq.length() - 1, seq.length());
            if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
                return "";
            } else {
                return decoySeq;
            }
        } else {
            char[] tempArray = seq.toCharArray();
            int idx = 0;
            while (idx < tempArray.length - 1) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            }
            String decoySeq = String.valueOf(tempArray);
            if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
                return "";
            } else {
                return decoySeq;
            }
        }
    }
}
