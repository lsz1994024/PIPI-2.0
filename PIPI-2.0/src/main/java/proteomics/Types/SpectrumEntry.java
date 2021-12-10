package proteomics.Types;

import java.util.*;

public class SpectrumEntry {
    public final int scanNum;
    public final double precursorMz;
    public final double precursorMass;
    public final int precursorCharge;
    public final TreeMap<Double, Double> plMap;
    private TreeMap<Double, Double> plMapWithVirtualPeaks = null;
    public final SparseVector plMapXcorr;
    private final String toString;

    public SpectrumEntry(int scanNum, double precursorMz, double precursorMass, int precursorCharge, TreeMap<Double, Double> plMap, SparseVector plMapXcorr) {
        this.scanNum = scanNum;
        this.precursorMz = precursorMz;
        this.precursorMass = precursorMass;
        this.precursorCharge = precursorCharge;
        this.plMap = plMap;
        this.plMapXcorr = plMapXcorr;
        toString = this.scanNum + " (charge = " + this.precursorCharge + ", mass = " + this.precursorMass + ", peak_num = " + this.plMap.size() + ")";
    }

    public String toString() {
        return toString;
    }

    public void setPlMapWithVirtualPeaks(TreeMap<Double, Double> input) {
        plMapWithVirtualPeaks = input;
    }

    public TreeMap<Double, Double> getPlMapWithVirtualPeaks() {
        return plMapWithVirtualPeaks;
    }
}