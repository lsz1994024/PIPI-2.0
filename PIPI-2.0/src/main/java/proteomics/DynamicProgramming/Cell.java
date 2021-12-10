package proteomics.DynamicProgramming;

import proteomics.Types.Coordinate;

public class Cell implements Comparable<Cell> {

    public final double v;
    public final Coordinate direction;

    public Cell(double v, Coordinate direction) {
        this.v = v;
        this.direction =  direction;
    }

    public int compareTo(Cell other) {
        if (this.v > other.v) {
            return 1;
        } else if (this.v < other.v) {
            return -1;
        } else {
            return 0;
        }
    }

    public String toString() {
        return String.format("%.1f", v) + "@" + direction;
    }
}