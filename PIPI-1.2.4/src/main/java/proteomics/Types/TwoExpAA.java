package proteomics.Types;


import java.util.Map;

public class TwoExpAA implements Comparable<TwoExpAA> {

    private final ExpAA[] twoExpAa;
    private String toString;
    private int hashCode;
    private final String aaString;
    private final float totalIntensity;
    private int regionIdx;
    private final float fidelityIndex;
    private Map<String, Float> massTable;

    public TwoExpAA(ExpAA aa1, ExpAA aa2, Map<String, Float> massTable) {
        twoExpAa = new ExpAA[]{aa1, aa2};
        toString = twoExpAa[0].toString() + "-" + twoExpAa[1].toString();
        hashCode = toString.hashCode();
        this.massTable = massTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : twoExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        float intensity = twoExpAa[0].getHeadIntensity();
        for (ExpAA aa : twoExpAa) {
            intensity += aa.getTailIntensity();
        }
//        totalIntensity = intensity;
        fidelityIndex = calFidelityIndex(twoExpAa);
        totalIntensity = intensity;
    }

    public float calFidelityIndex(ExpAA[] twoExpAa) {
        Map<String, Float> massTable = this.massTable;

        float peak1 = twoExpAa[0].getHeadLocation();
        float peak2 = twoExpAa[1].getHeadLocation();
        float peak3 = twoExpAa[1].getTailLocation();
        float mass1 = massTable.get(String.valueOf(twoExpAa[0].getAA()));
        float mass2 = massTable.get(String.valueOf(twoExpAa[1].getAA()));

        // first peak is true
        float theoPeak12 = peak1 + mass1;
        float theoPeak13 = peak1 + mass1 + mass2;
        float diff12 = peak2 - theoPeak12;
        float diff13 = peak3 - theoPeak13;
        float value12 = calGaussian(2, diff12*100);
        float value13 = calGaussian(4, diff13*100);

        //second peak is true
        float theoPeak21 = peak2 - mass1;
        float theoPeak23 = peak2 + mass2;
        float diff21 = peak1 - theoPeak21;
        float diff23 = peak3 - theoPeak23;
        float value21 = calGaussian(2, diff21*100);
        float value23 = calGaussian(2, diff23*100);

        //third peak is true
        float theoPeak31 = peak3 - mass2 - mass1;
        float theoPeak32 = peak3 - mass2;
        float diff31 = peak1 - theoPeak31;
        float diff32 = peak2 - theoPeak32;
        float value31 = calGaussian(4, diff31*100);
        float value32 = calGaussian(2, diff32*100);

        float value = (value12+value13+value21+value23+value31+value32)/6;

        return value;
    }

    public float calGaussian (double theta, double x) {
        float nominator = (float) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        float denominator = (float) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
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

    public float getTotalIntensity() {
        return totalIntensity;
    }

    public float getHeadLocation() {
        return twoExpAa[0].getHeadLocation();
    }

    public float getFidelityIndex() { return fidelityIndex; }

    public float getTailLocation() {
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
