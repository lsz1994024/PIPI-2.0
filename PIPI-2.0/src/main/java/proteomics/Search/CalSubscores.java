package proteomics.Search;


import proteomics.Types.FinalResultEntry;
import proteomics.Types.Peptide;
import proteomics.Types.SpectrumEntry;

import java.util.*;

public class CalSubscores {

    public CalSubscores(List<FinalResultEntry> psmList, Map<Integer, SpectrumEntry> numSpectrumMap, double ms2Tolerance) {
        for (FinalResultEntry psm : psmList) {
            SpectrumEntry spectrum = numSpectrumMap.get(psm.getScanNum());
            TreeMap<Double, Double> expPl = spectrum.plMap;
            Peptide peptide = psm.getPeptide();
            double[][] ionMatrix = peptide.getIonMatrix();
            int precursorCharge = spectrum.precursorCharge;

            double matchedPeakIntensity = 0;
            int maxRow = Math.min(ionMatrix.length, 2 * (precursorCharge - 1));
            if (precursorCharge == 1) {
                maxRow = 2;
            }
            int totalIonNum = ionMatrix[0].length * maxRow;
            Double[] intensityArray = expPl.values().toArray(new Double[expPl.size()]);
            Arrays.sort(intensityArray, Collections.reverseOrder());
            double intensityT = 0;
            if (totalIonNum < intensityArray.length) {
                intensityT = intensityArray[totalIonNum];
            }
            int matchedHighestPeakNum = 0;
            double totalIntensity = 0;
            for (int i = 0; i < maxRow; ++i) {
                for (int j = 0; j < ionMatrix[0].length; ++j) {
                    for (double mz : expPl.keySet()) {
                        double intensity = expPl.get(mz);
                        totalIntensity += intensity;
                        if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                            if (intensity > intensityT) {
                                ++matchedHighestPeakNum;
                            }
                            matchedPeakIntensity += intensity;
                            break;
                        }
                    }
                }
            }

            psm.setIonFrac(matchedPeakIntensity / totalIntensity);
            psm.setMatchedHighestIntensityFrac((double) matchedHighestPeakNum / (double) totalIonNum);
        }
    }
}
