package qut.edu.au.pqcsebase.projects.IBEKS;


public class Params {

    public static final int flag = 1; // 0 for trapGenFirst, 1 for second

    public static final long n = 5; // Dimension of the lattice
    public static final long q = 17; // A prime modulus

    public static final long N = 3; //keywords set n > N + 1

    public static final long n_N = n - N - 1; // n - N - 1

    public static final double sigma = 1000000.0; // Gaussian parameter

    public static final String[] keywords = new String[]{"doctor", "engineer", "teacher", "lawyer", "nurse"};
}