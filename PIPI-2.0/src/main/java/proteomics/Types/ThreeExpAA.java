package proteomics.Types;


import java.util.Map;

public class ThreeExpAA implements Comparable<ThreeExpAA> {

    private final ExpAA[] threeExpAa;
    private String toString;
    private int hashCode;
    private final String aaString;
    private final double totalIntensity;
    private int regionIdx;
    private final double fidelityIndex;
    private Map<Character, Double> massTable;

    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, Map<Character, Double> massTable) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3};
        toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();
        this.massTable = massTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }

        fidelityIndex = calFidelityIndex(threeExpAa);
        totalIntensity = intensity;
    }

    public double calFidelityIndex(ExpAA[] threeExpAa) {
        Map<Character, Double> massTable = this.massTable;

        double peak1 = threeExpAa[0].getHeadLocation();
        double peak2 = threeExpAa[1].getHeadLocation();
        double peak3 = threeExpAa[2].getHeadLocation();
        double peak4 = threeExpAa[2].getTailLocation();
        double mass1 = massTable.get(String.valueOf(threeExpAa[0].getAA()));
        double mass2 = massTable.get(String.valueOf(threeExpAa[1].getAA()));
        double mass3 = massTable.get(String.valueOf(threeExpAa[2].getAA()));

        // first peak is true
        double theoPeak12 = peak1 + mass1;
        double theoPeak13 = peak1 + mass1 + mass2;
        double theoPeak14 = peak1 + mass1 + mass2 + mass3;
        double diff12 = peak2 - theoPeak12;
        double diff13 = peak3 - theoPeak13;
        double diff14 = peak4 - theoPeak14;
        double value12 = calGaussian(2, diff12*100);
        double value13 = calGaussian(4, diff13*100);
        double value14 = calGaussian(6, diff14*100);

        //second peak is true
        double theoPeak21 = peak2 - mass1;
        double theoPeak23 = peak2 + mass2;
        double theoPeak24 = peak2 + mass2 + mass3;
        double diff21 = peak1 - theoPeak21;
        double diff23 = peak3 - theoPeak23;
        double diff24 = peak4 - theoPeak24;
        double value21 = calGaussian(2, diff21*100);
        double value23 = calGaussian(2, diff23*100);
        double value24 = calGaussian(4, diff24*100);

        //third peak is true
        double theoPeak31 = peak3 - mass2 - mass1;
        double theoPeak32 = peak3 - mass2;
        double theoPeak34 = peak3 + mass3;
        double diff31 = peak1 - theoPeak31;
        double diff32 = peak2 - theoPeak32;
        double diff34 = peak4 - theoPeak34;
        double value31 = calGaussian(4, diff31*100);
        double value32 = calGaussian(2, diff32*100);
        double value34 = calGaussian(2, diff34*100);

        //fourth peak is true
        double theoPeak41 = peak4 - mass3 - mass2 - mass1;
        double theoPeak42 = peak4 - mass3 - mass2;
        double theoPeak43 = peak4 - mass3;
        double diff41 = peak1 - theoPeak41;
        double diff42 = peak2 - theoPeak42;
        double diff43 = peak3 - theoPeak43;
        double value41 = calGaussian(6, diff41*100);
        double value42 = calGaussian(4, diff42*100);
        double value43 = calGaussian(2, diff43*100);


        double value = (value12+value13+value14+value21+value23+value24+value31+value32+value34+value41+value42+value43)/12;

        return value;
    }

    public double calGaussian (double theta, double x) {
        double nominator = (double) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        double denominator = (double) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
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

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getFidelityIndex() { return fidelityIndex; }

    public double getHeadLocation() {
        return threeExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
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
