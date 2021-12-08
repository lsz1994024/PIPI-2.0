package proteomics.Types;


import java.util.Map;

public class ThreeExpAA implements Comparable<ThreeExpAA> {

    private final ExpAA[] threeExpAa;
    private String toString;
    private int hashCode;
    private final String aaString;
    private final float totalIntensity;
    private int regionIdx;
    private final float fidelityIndex;
    private Map<String, Float> massTable;

    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, Map<String, Float> massTable) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3};
        toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();
        this.massTable = massTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        float intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }

        fidelityIndex = calFidelityIndex(threeExpAa);
        totalIntensity = intensity;
    }

    public float calFidelityIndex(ExpAA[] threeExpAa) {
        Map<String, Float> massTable = this.massTable;

        float peak1 = threeExpAa[0].getHeadLocation();
        float peak2 = threeExpAa[1].getHeadLocation();
        float peak3 = threeExpAa[2].getHeadLocation();
        float peak4 = threeExpAa[2].getTailLocation();
        float mass1 = massTable.get(String.valueOf(threeExpAa[0].getAA()));
        float mass2 = massTable.get(String.valueOf(threeExpAa[1].getAA()));
        float mass3 = massTable.get(String.valueOf(threeExpAa[2].getAA()));

        // first peak is true
        float theoPeak12 = peak1 + mass1;
        float theoPeak13 = peak1 + mass1 + mass2;
        float theoPeak14 = peak1 + mass1 + mass2 + mass3;
        float diff12 = peak2 - theoPeak12;
        float diff13 = peak3 - theoPeak13;
        float diff14 = peak4 - theoPeak14;
        float value12 = calGaussian(2, diff12*100);
        float value13 = calGaussian(4, diff13*100);
        float value14 = calGaussian(6, diff14*100);

        //second peak is true
        float theoPeak21 = peak2 - mass1;
        float theoPeak23 = peak2 + mass2;
        float theoPeak24 = peak2 + mass2 + mass3;
        float diff21 = peak1 - theoPeak21;
        float diff23 = peak3 - theoPeak23;
        float diff24 = peak4 - theoPeak24;
        float value21 = calGaussian(2, diff21*100);
        float value23 = calGaussian(2, diff23*100);
        float value24 = calGaussian(4, diff24*100);

        //third peak is true
        float theoPeak31 = peak3 - mass2 - mass1;
        float theoPeak32 = peak3 - mass2;
        float theoPeak34 = peak3 + mass3;
        float diff31 = peak1 - theoPeak31;
        float diff32 = peak2 - theoPeak32;
        float diff34 = peak4 - theoPeak34;
        float value31 = calGaussian(4, diff31*100);
        float value32 = calGaussian(2, diff32*100);
        float value34 = calGaussian(2, diff34*100);

        //fourth peak is true
        float theoPeak41 = peak4 - mass3 - mass2 - mass1;
        float theoPeak42 = peak4 - mass3 - mass2;
        float theoPeak43 = peak4 - mass3;
        float diff41 = peak1 - theoPeak41;
        float diff42 = peak2 - theoPeak42;
        float diff43 = peak3 - theoPeak43;
        float value41 = calGaussian(6, diff41*100);
        float value42 = calGaussian(4, diff42*100);
        float value43 = calGaussian(2, diff43*100);


        float value = (value12+value13+value14+value21+value23+value24+value31+value32+value34+value41+value42+value43)/12;

        return value;
    }

    public float calGaussian (double theta, double x) {
        float nominator = (float) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        float denominator = (float) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        return nominator/denominator;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ThreeExpAA) && (this.hashCode() == other.hashCode());
    }

    public void setTheoLocation(int i, int theoLoc) {
        threeExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int compareTo(ThreeExpAA other) {
        if (this.threeExpAa[0].getHeadLocation() > other.getExpAAs()[0].getHeadLocation()) {
            return 1;
        } else if (this.threeExpAa[0].getHeadLocation() < other.getExpAAs()[0].getHeadLocation()) {
            return -1;
        } else {
            return 0;
        }
    }

    public ExpAA[] getExpAAs() {
        return threeExpAa;
    }

    public String getAAString() {
        return aaString;
    }

    public float getTotalIntensity() {
        return totalIntensity;
    }

    public float getFidelityIndex() { return fidelityIndex; }

    public float getHeadLocation() {
        return threeExpAa[0].getHeadLocation();
    }

    public float getTailLocation() {
        return threeExpAa[threeExpAa.length - 1].getTailLocation();
    }

    public ThreeExpAA clone() {
        return new ThreeExpAA(threeExpAa[0].clone(), threeExpAa[1].clone(), threeExpAa[2].clone(), massTable);
    }

    public int size() {
        return threeExpAa.length;
    }

    public ExpAA get(int i) {
        return threeExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
