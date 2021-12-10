package proteomics.Spectrum;

import proteomics.Library.MassTool;
import proteomics.Types.SparseVector;

import java.util.*;

public class PreSpectrum {

    private static final double defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final double doubleZero = 1e-6f;
    private static final int xcprrOffset = 75;

    private final MassTool massToolObj;
    private final Map<Character, Double> massTable;

    public PreSpectrum(MassTool massToolObj) {
        this.massToolObj = massToolObj;
        massTable = massToolObj.getMassTable();
    }

    public TreeMap<Double, Double> preSpectrum (Map<Double, Double> peaksMap, double precursorMass, int precursorCharge, double ms2Tolerance) {
        // remove precursor peak from spectrum
        TreeMap<Double, Double> temp = removePrecursorPeak(peaksMap, precursorMass, precursorCharge, ms2Tolerance);

        // reduce noise
        TreeMap<Double, Double> deionisedPlMap = deNoise(new TreeMap<>(temp.subMap(0.0, precursorMass)));

        // normalize
        return normalizeSpec(deionisedPlMap);
    }

    public SparseVector prepareXcorr(TreeMap<Double, Double> plMap, double precursorMass) {
        double[] plArray = digitizeSpec(plMap, precursorMass);

        SparseVector xcorrPl = new SparseVector();
        double mySum = 0;
        int offsetRange = 2 * xcprrOffset + 1;
        for (int i = 0; i < xcprrOffset; ++i) {
            mySum += plArray[i];
        }

        double factor = 1 / (double) (offsetRange - 1);
        for (int i = xcprrOffset; i < plArray.length + xcprrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            double temp = (plArray[i - xcprrOffset] - (mySum - plArray[i - xcprrOffset]) * factor);
            if (Math.abs(temp) > doubleZero) {
                xcorrPl.put(i - xcprrOffset, temp);
            }
        }

        return xcorrPl;
    }

    private TreeMap<Double, Double> removePrecursorPeak(Map<Double, Double> peakMap, double precursorMass, int precursorCharge, double ms2Tolerance) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();

        for (double mz : peakMap.keySet()) {
            for (int charge = precursorCharge; charge > 0; --charge) {
                double temp = (precursorMass + charge * massTable.get("PROTON")) / charge;
                if ((peakMap.get(mz) > doubleZero) && (Math.abs(peakMap.get(mz) - temp) > ms2Tolerance)) {
                    mzIntensityMap.put((double) mz, peakMap.get(mz).doubleValue());
                }
            }
        }

        return mzIntensityMap;
    }

    private TreeMap<Double, Double> deNoise(TreeMap<Double, Double> plMap) {
        // denoise
        TreeMap<Double, Double> denoisedPlMap = new TreeMap<>();
        double minMz = plMap.firstKey();
        double maxMz = plMap.lastKey();
        double windowSize = (plMap.lastKey() - plMap.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            double leftMz = Math.min(minMz + i * windowSize, maxMz);
            double rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Double, Double> subPlMap;
            if (rightMz < maxMz) {
                subPlMap = plMap.subMap(leftMz, true, rightMz, false);
            } else {
                subPlMap = plMap.subMap(leftMz, true, rightMz, true);
            }

            if (subPlMap.size() > 9) {
                double noiseIntensity = estimateNoiseIntensity(subPlMap);
                for (double mz : subPlMap.keySet()) {
                    if (subPlMap.get(mz) > noiseIntensity) {
                        denoisedPlMap.put(mz, subPlMap.get(mz));
                    }
                }
            } else {
                for (double mz : subPlMap.keySet()) {
                    denoisedPlMap.put(mz, subPlMap.get(mz));
                }
            }
        }

        return denoisedPlMap;
    }

    private double estimateNoiseIntensity(Map<Double, Double> pl) {
        Set<Double> intensitySet = new HashSet<>();
        for (double intensity : pl.values()) {
            intensitySet.add(intensity);
        }
        Double[] uniqueIntensityVector = intensitySet.toArray(new Double[intensitySet.size()]);
        Arrays.sort(uniqueIntensityVector);
        double[] cum = new double[uniqueIntensityVector.length];
        for (int i = 0; i < uniqueIntensityVector.length; ++i) {
            for (double intensity : pl.values()) {
                if (intensity <= uniqueIntensityVector[i]) {
                    ++cum[i];
                }
            }
        }
        double[][] diff = new double[2][uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = uniqueIntensityVector[i + 1] - uniqueIntensityVector[i];
        }
        double[] diff2 = new double[uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff2[i] = diff[0][i] / (diff[1][i] + doubleZero);
        }
        double maxValue = 0;
        int maxIdx = 0;
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            if (diff2[i] > maxValue) {
                maxValue = diff2[i];
                maxIdx = i;
            }
        }

        return uniqueIntensityVector[maxIdx];
    }

    private TreeMap<Double, Double> normalizeSpec(TreeMap<Double, Double> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            if (plMap.get(mz) > doubleZero) {
                double sqrtIntensity = (double) Math.sqrt(plMap.get(mz));
                sqrtPlMap.put(mz, sqrtIntensity);
            }
        }

        // divide the spectrum into 10 windows and normalize each windows to defaultIntensity
        TreeMap<Double, Double> windowedPlMap = new TreeMap<>();
        double minMz = sqrtPlMap.firstKey();
        double maxMz = sqrtPlMap.lastKey();
        double windowSize = (maxMz - minMz) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            double leftMz = Math.min(minMz + i * windowSize, maxMz);
            double rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Double, Double> subMap;
            if (rightMz < maxMz) {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, false);
            } else {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, true);
            }
            if (!subMap.isEmpty()) {
                Double[] intensityArray = subMap.values().toArray(new Double[subMap.size()]);
                Arrays.sort(intensityArray);
                double temp1 = defaultIntensity / intensityArray[intensityArray.length - 1];
                double temp2 = (double) 0.05 * intensityArray[intensityArray.length - 1];
                for (double mz : subMap.keySet()) {
                    if (subMap.get(mz) > temp2) {
                        windowedPlMap.put(mz, subMap.get(mz) * temp1);
                    }
                }
            }
        }

        return windowedPlMap;
    }

    private double[] digitizeSpec(TreeMap<Double, Double> pl, double precursorMass) {
        double[] digitizedPl = new double[massToolObj.mzToBin(precursorMass) + 1];
        Set<Double> mzSet = pl.keySet();
        for (double mz : mzSet) {
            int idx = massToolObj.mzToBin(mz);
            digitizedPl[idx] = Math.max(pl.get(mz), digitizedPl[idx]);
        }

        return digitizedPl;
    }
}
