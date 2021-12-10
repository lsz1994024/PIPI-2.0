package proteomics.Types;

public class ExpAA implements Comparable<ExpAA> {

    private final char aa;
    private final double headLocation;
    private final double tailLocation;
    private final double headIntensity;
    private final double tailIntensity;
    private final double totalHalfIntensity;
    private int theoLocation; // starts from 0
    private String toString;
    private int hashCode;

    public ExpAA(char aa, double headLocation, double tailLocation, double headIntensity, double tailIntensity, int theoLocation) {
        this.aa = aa;
        this.headLocation = headLocation;
        this.tailLocation = tailLocation;
        this.headIntensity = headIntensity;
        this.tailIntensity = tailIntensity;
        this.totalHalfIntensity = (headIntensity + tailIntensity) / 2;
        this.theoLocation = -1;
        toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public char getAA() {
        return aa;
    }

    protected void setTheoLocation(int theo) {
        theoLocation = theo;
        // update toString and hashCode
        toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public int getTheoLocation() {
        return theoLocation;
    }

    public int compareTo(ExpAA other) {
        if (headLocation > other.getHeadLocation()) {
            return 1;
        } else if (headLocation < other.getHeadLocation()) {
            return -1;
        } else {
            return 0;
        }
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ExpAA) {
            ExpAA temp = (ExpAA) other;
            return this.hashCode() == temp.hashCode();
        } else {
            return false;
        }
    }

    public boolean approximateEquals(ExpAA other, double tolerance) {
        return ((this.aa == other.aa) && (this.theoLocation == other.theoLocation) && (Math.abs(this.headLocation - other.headLocation) <= tolerance));
    }

    public ExpAA clone() {
        return new ExpAA(aa, headLocation, tailLocation, headIntensity, tailIntensity, theoLocation);
    }

    public double getHeadLocation() {
        return headLocation;
    }

    public double getTailLocation() {
        return tailLocation;
    }

    public double getHeadIntensity() {
        return headIntensity;
    }

    public double getTailIntensity() {
        return tailIntensity;
    }

    public double getTotalHalfIntensity() {
        return totalHalfIntensity;
    }
}