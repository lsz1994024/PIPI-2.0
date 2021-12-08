package proteomics.Types;

//import ProteomicsLibrary.MassTool;

import java.util.Map;

public class FourExpAA implements Comparable<FourExpAA> {

    private final ExpAA[] fourExpAa;
    private int hashCode;
    private String toString;
    private final String aaString;
    private final float totalIntensity;
    private final float fidelityIndex;
    private int regionIdx;
    private Map<String, Float> massTable;

    public FourExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, ExpAA aa4, Map<String, Float> MassTable) {
        fourExpAa = new ExpAA[]{aa1, aa2, aa3, aa4};
        toString = fourExpAa[0].toString() + "-" + fourExpAa[1].toString() + "-" + fourExpAa[2].toString() + "-" + fourExpAa[3].toString();
        hashCode = toString.hashCode();
        this.massTable = MassTable;

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : fourExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        float intensity = fourExpAa[0].getHeadIntensity();
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

    public float calFidelityIndex(ExpAA[] fourExpAa) {
        Map<String, Float> massTable = this.massTable;

        float peak1 = fourExpAa[0].getHeadLocation();
        float peak2 = fourExpAa[1].getHeadLocation();
        float peak3 = fourExpAa[2].getHeadLocation();
        float peak4 = fourExpAa[3].getHeadLocation();
        float peak5 = fourExpAa[3].getTailLocation();
        float mass1 = massTable.get(String.valueOf(fourExpAa[0].getAA()));
        float mass2 = massTable.get(String.valueOf(fourExpAa[1].getAA()));
        float mass3 = massTable.get(String.valueOf(fourExpAa[2].getAA()));
        float mass4 = massTable.get(String.valueOf(fourExpAa[3].getAA()));

        // first peak is true
        float theoPeak12 = peak1 + mass1;
        float theoPeak13 = peak1 + mass1 + mass2;
        float theoPeak14 = peak1 + mass1 + mass2 + mass3;
        float theoPeak15 = peak1 + mass1 + mass2 + mass3 + mass4;
        float diff12 = peak2 - theoPeak12;
        float diff13 = peak3 - theoPeak13;
        float diff14 = peak4 - theoPeak14;
        float diff15 = peak5 - theoPeak15;
        float value12 = calGaussian(2, diff12*100);
        float value13 = calGaussian(4, diff13*100);
        float value14 = calGaussian(6, diff14*100);
        float value15 = calGaussian(8, diff15*100);

        //second peak is true
        float theoPeak21 = peak2 - mass1;
        float theoPeak23 = peak2 + mass2;
        float theoPeak24 = peak2 + mass2 + mass3;
        float theoPeak25 = peak2 + mass2 + mass3 + mass4;
        float diff21 = peak1 - theoPeak21;
        float diff23 = peak3 - theoPeak23;
        float diff24 = peak4 - theoPeak24;
        float diff25 = peak5 - theoPeak25;
        float value21 = calGaussian(2, diff21*100);
        float value23 = calGaussian(2, diff23*100);
        float value24 = calGaussian(4, diff24*100);
        float value25 = calGaussian(6, diff25*100);

        //third peak is true
        float theoPeak31 = peak3 - mass2 - mass1;
        float theoPeak32 = peak3 - mass2;
        float theoPeak34 = peak3 + mass3;
        float theoPeak35 = peak3 + mass3 + mass4;
        float diff31 = peak1 - theoPeak31;
        float diff32 = peak2 - theoPeak32;
        float diff34 = peak4 - theoPeak34;
        float diff35 = peak5 - theoPeak35;
        float value31 = calGaussian(4, diff31*100);
        float value32 = calGaussian(2, diff32*100);
        float value34 = calGaussian(2, diff34*100);
        float value35 = calGaussian(4, diff35*100);

        //fourth peak is true
        float theoPeak41 = peak4 - mass3 - mass2 - mass1;
        float theoPeak42 = peak4 - mass3 - mass2;
        float theoPeak43 = peak4 - mass3;
        float theoPeak45 = peak4 + mass4;
        float diff41 = peak1 - theoPeak41;
        float diff42 = peak2 - theoPeak42;
        float diff43 = peak3 - theoPeak43;
        float diff45 = peak5 - theoPeak45;
        float value41 = calGaussian(6, diff41*100);
        float value42 = calGaussian(4, diff42*100);
        float value43 = calGaussian(2, diff43*100);
        float value45 = calGaussian(2, diff45*100);

        //fifth peak is true
        float theoPeak51 = peak5 - mass4 - mass3 - mass2 - mass1;
        float theoPeak52 = peak5 - mass4 - mass3 - mass2;
        float theoPeak53 = peak5 - mass4 - mass3;
        float theoPeak54 = peak5 - mass4;
        float diff51 = peak1 - theoPeak51;
        float diff52 = peak2 - theoPeak52;
        float diff53 = peak3 - theoPeak53;
        float diff54 = peak5 - theoPeak54;
        float value51 = calGaussian(6, diff51*100);
        float value52 = calGaussian(4, diff52*100);
        float value53 = calGaussian(2, diff53*100);
        float value54 = calGaussian(2, diff54*100);

        float value = (value12+value13+value14+value15+value21+value23+value24+value25+value31+value32+value34+value35+value41+value42+value43+value45+value51+value52+value53+value54)/20;

        return value;
    }

    public float calGaussian (double theta, double x) {
        float nominator = (float) (Math.pow(Math.E, -(x*x)/(2*theta))/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
        float denominator = (float) (Math.pow(Math.E, 0)/(Math.sqrt(theta)*Math.sqrt(2*Math.PI)));
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

    public float getTotalIntensity() {
        return totalIntensity;
    }

    public float getFidelityIndex() { return fidelityIndex; }

    public float getHeadLocation() {
        return fourExpAa[0].getHeadLocation();
    }

    public float getTailLocation() {
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
