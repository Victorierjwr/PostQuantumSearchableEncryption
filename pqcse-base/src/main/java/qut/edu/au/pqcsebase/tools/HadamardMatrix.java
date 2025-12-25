package qut.edu.au.pqcsebase.tools;

import java.util.*;

public class HadamardMatrix {
    // Generate Hadamard matrix H_n (n must be power of two), entries are +1/-1
    public static int[][] hadamard(int n) {
        if (n <= 0 || (n & (n - 1)) != 0) {
            throw new IllegalArgumentException("n must be a power of two");
        }
        int[][] H = new int[n][n];
        H[0][0] = 1;
        for (int size = 1; size < n; size <<= 1) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    int v = H[i][j];
                    H[i][j + size] = v;      // top-right
                    H[i + size][j] = v;      // bottom-left
                    H[i + size][j + size] = -v; // bottom-right
                }
            }
        }
        return H;
    }

    public static void main(String[] args) {
        int n = 16;
        int[][] H = hadamard(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(H[i][j] + (j + 1 == n ? "" : " "));
            }
            System.out.println();
        }
    }
}
