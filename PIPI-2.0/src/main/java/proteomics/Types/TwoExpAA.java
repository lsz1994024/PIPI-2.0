package proteomics.Types;


import java.util.Map;

public class TwoExpAA implements Comparable<TwoExpAA> {

    private final ExpAA[] twoExpAa;
    private String toString;
    private int hashCode;
    private final String aaString;
    private final double
            totalIntensity;
    private int regionIdx;
    private final double
            fidelityIndex;
    private Map<Character, Double> massTable;

    public TwoExpAA(ExpAA aa1, ExpAA aa2, Map<Character, Double> massTable) {
        twoExpAa = new ExpAA[]{aa1, aa2};
        toString = twoExpAa[0].toString() + "-" + twoExpAa[1].toString();
        hashCode = toString.hashCode();
        this.massTable = massTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : twoExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        double
                intensity = twoExpAa[0].getHeadIntensity();
        for (ExpAA aa : twoExpAa) {
            intensity += aa.getTailIntensity();
        }
//        totalIntensity = intensity;
        fidelityIndex = calFidelityIndex(twoExpAa);
        totalIntensity = intensity;
    }

    public double
    calFidelityIndex(ExpAA[] twoExpAa) {
        Map<Character, Double> massTable = this.massTable;

        double
                peak1 = twoExpAa[0].getHeadLocation();
        double
                peak2 = twoExpAa[1].getHeadLocation();
        double
                peak3 = twoExpAa[1].getTailLocation();
        double
                mass1 = massTable.get(String.valueOf(twoExpAa[0].getAA()));
        double
                mass2 = massTable.get(String.valueOf(twoExpAa[1].getAA()));

        // first peak is true
        double
                theoPeak12 = peak1 + mass1;
        double
                theoPeak13 = peak1 + mass1 + mass2;
        double
                diff12 = peak2 - theoPeak12;
        double
                diff13 = peak3 - theoPeak13;
        double
                value12 = calGaussian(2, diff12*100);
        double
                value13 = calGaussian(4, diff13*100);

        //second peak is true
        double
                theoPeak21 = peak2 - mass1;
        double
                theoPeak23 = peak2 + mass2;
        double
                diff21 = peak1 - theoPeak21;
        double
                diff23 = peak3 - theoPeak23;
        double
                value21 = calGaussian(2, diff21*100);
        double
                value23 = calGaussian(2, diff23*100);

        //third peak is true
        double
                theoPeak31 = peak3 - mass2 - mass1;
        double
                theoPeak32 = peak3 - mass2;
        double
                diff31 = peak1 - theoPeak31;
        double
                diff32 = peak2 - theoPeak32;
        double
                value31 = calGaussian(4, diff31*100);
        double
                value32 = calGaussian(2, diff32*100);

        double
                value = (value12+value13+value21+value23+value31+value32)/6;

        return value;
    }

    public double
    calGaussian (double theta, double x) {
        double
                nominator = (double
                ) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        double
                denominator = (double
                ) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        return nominator/denominator;
    }

    public boolean equals(Object other) {
        return (other instanceof TwoExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(TwoExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public int hashCode() {
        return hashCode;
    }

    public void setTheoLocation(int i, int theoLoc) {
        twoExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        String toString = twoExpAa[0].toString() + "-" + twoExpAa[1].toString() + "-" + twoExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    public int compareTo(TwoExpAA other) {
        return Double.compare(twoExpAa[0].getHeadLocation(), other.twoExpAa[0].getHeadLocation());
    }

    public ExpAA[] getExpAAs() {
        return twoExpAa;
    }

    public String getAAString() {
        return aaString;
    }

    public double
    getTotalIntensity() {
        return totalIntensity;
    }

    public double
    getHeadLocation() {
        return twoExpAa[0].getHeadLocation();
    }

    public double
    getFidelityIndex() { return fidelityIndex; }

    public double
    getTailLocation() {
        return twoExpAa[twoExpAa.length - 1].getTailLocation();
    }

    public TwoExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new TwoExpAA(twoExpAa[0].clone(), twoExpAa[1].clone(), massTable);
    }

    public int size() {
        return twoExpAa.length;
    }

    public ExpAA get(int i) {
        return twoExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
