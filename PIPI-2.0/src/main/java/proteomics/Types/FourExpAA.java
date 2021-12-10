package proteomics.Types;

//import ProteomicsLibrary.MassTool;

import java.util.Map;

public class FourExpAA implements Comparable<FourExpAA> {

    private final ExpAA[] fourExpAa;
    private int hashCode;
    private String toString;
    private final String aaString;
    private final double totalIntensity;
    private final double fidelityIndex;
    private int regionIdx;
    private Map<Character, Double> massTable;

    public FourExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, ExpAA aa4, Map<Character, Double> MassTable) {
        fourExpAa = new ExpAA[]{aa1, aa2, aa3, aa4};
        toString = fourExpAa[0].toString() + "-" + fourExpAa[1].toString() + "-" + fourExpAa[2].toString() + "-" + fourExpAa[3].toString();
        hashCode = toString.hashCode();
        this.massTable = MassTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : fourExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        double intensity = fourExpAa[0].getHeadIntensity();
        for (ExpAA aa : fourExpAa) {
            intensity += aa.getTailIntensity();
        }
//        totalIntensity = intensity;
        fidelityIndex = calFidelityIndex(fourExpAa);
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public double calFidelityIndex(ExpAA[] fourExpAa) {
        Map<Character, Double> massTable = this.massTable;

        double peak1 = fourExpAa[0].getHeadLocation();
        double peak2 = fourExpAa[1].getHeadLocation();
        double peak3 = fourExpAa[2].getHeadLocation();
        double peak4 = fourExpAa[3].getHeadLocation();
        double peak5 = fourExpAa[3].getTailLocation();
        double mass1 = massTable.get(String.valueOf(fourExpAa[0].getAA()));
        double mass2 = massTable.get(String.valueOf(fourExpAa[1].getAA()));
        double mass3 = massTable.get(String.valueOf(fourExpAa[2].getAA()));
        double mass4 = massTable.get(String.valueOf(fourExpAa[3].getAA()));

        // first peak is true
        double theoPeak12 = peak1 + mass1;
        double theoPeak13 = peak1 + mass1 + mass2;
        double theoPeak14 = peak1 + mass1 + mass2 + mass3;
        double theoPeak15 = peak1 + mass1 + mass2 + mass3 + mass4;
        double diff12 = peak2 - theoPeak12;
        double diff13 = peak3 - theoPeak13;
        double diff14 = peak4 - theoPeak14;
        double diff15 = peak5 - theoPeak15;
        double value12 = calGaussian(2, diff12*100);
        double value13 = calGaussian(4, diff13*100);
        double value14 = calGaussian(6, diff14*100);
        double value15 = calGaussian(8, diff15*100);

        //second peak is true
        double theoPeak21 = peak2 - mass1;
        double theoPeak23 = peak2 + mass2;
        double theoPeak24 = peak2 + mass2 + mass3;
        double theoPeak25 = peak2 + mass2 + mass3 + mass4;
        double diff21 = peak1 - theoPeak21;
        double diff23 = peak3 - theoPeak23;
        double diff24 = peak4 - theoPeak24;
        double diff25 = peak5 - theoPeak25;
        double value21 = calGaussian(2, diff21*100);
        double value23 = calGaussian(2, diff23*100);
        double value24 = calGaussian(4, diff24*100);
        double value25 = calGaussian(6, diff25*100);

        //third peak is true
        double theoPeak31 = peak3 - mass2 - mass1;
        double theoPeak32 = peak3 - mass2;
        double theoPeak34 = peak3 + mass3;
        double theoPeak35 = peak3 + mass3 + mass4;
        double diff31 = peak1 - theoPeak31;
        double diff32 = peak2 - theoPeak32;
        double diff34 = peak4 - theoPeak34;
        double diff35 = peak5 - theoPeak35;
        double value31 = calGaussian(4, diff31*100);
        double value32 = calGaussian(2, diff32*100);
        double value34 = calGaussian(2, diff34*100);
        double value35 = calGaussian(4, diff35*100);

        //fourth peak is true
        double theoPeak41 = peak4 - mass3 - mass2 - mass1;
        double theoPeak42 = peak4 - mass3 - mass2;
        double theoPeak43 = peak4 - mass3;
        double theoPeak45 = peak4 + mass4;
        double diff41 = peak1 - theoPeak41;
        double diff42 = peak2 - theoPeak42;
        double diff43 = peak3 - theoPeak43;
        double diff45 = peak5 - theoPeak45;
        double value41 = calGaussian(6, diff41*100);
        double value42 = calGaussian(4, diff42*100);
        double value43 = calGaussian(2, diff43*100);
        double value45 = calGaussian(2, diff45*100);

        //fifth peak is true
        double theoPeak51 = peak5 - mass4 - mass3 - mass2 - mass1;
        double theoPeak52 = peak5 - mass4 - mass3 - mass2;
        double theoPeak53 = peak5 - mass4 - mass3;
        double theoPeak54 = peak5 - mass4;
        double diff51 = peak1 - theoPeak51;
        double diff52 = peak2 - theoPeak52;
        double diff53 = peak3 - theoPeak53;
        double diff54 = peak5 - theoPeak54;
        double value51 = calGaussian(6, diff51*100);
        double value52 = calGaussian(4, diff52*100);
        double value53 = calGaussian(2, diff53*100);
        double value54 = calGaussian(2, diff54*100);

        double value = (value12+value13+value14+value15+value21+value23+value24+value25+value31+value32+value34+value35+value41+value42+value43+value45+value51+value52+value53+value54)/20;

        return value;
    }

    public double calGaussian (double theta, double x) {
        double nominator = (double) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        double denominator = (double) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        return nominator/denominator;
    }

    public boolean equals(Object other) {
        return (other instanceof FourExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(FourExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public void setTheoLocation(int i, int theoLoc) {
        fourExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        String toString = fourExpAa[0].toString() + "-" + fourExpAa[1].toString() + "-" + fourExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    public int compareTo(FourExpAA other) {
        return Double.compare(fourExpAa[0].getHeadLocation(), other.fourExpAa[0].getHeadLocation());
    }

    public ExpAA[] getExpAAs() {
        return fourExpAa;
    }

    public String getAaString() {
        return aaString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getFidelityIndex() { return fidelityIndex; }

    public double getHeadLocation() {
        return fourExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
        return fourExpAa[fourExpAa.length - 1].getTailLocation();
    }

    public FourExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new FourExpAA(fourExpAa[0].clone(), fourExpAa[1].clone(), fourExpAa[2].clone(), fourExpAa[3].clone(), massTable);
    }

    public int size() {
        return fourExpAa.length;
    }

    public ExpAA get(int i) {
        return fourExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
