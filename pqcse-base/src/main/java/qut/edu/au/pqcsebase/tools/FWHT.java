package qut.edu.au.pqcsebase.tools;

import java.util.*;

public class FWHT {
    // XOR-FWHT in-place. If inverse=true, does inverse transform (divide by n at end).
    // Requires n = power of two.
    public static void fwhtXor(long[] a, boolean inverse) {
        int n = a.length;
        for (int len = 1; len < n; len <<= 1) {
            for (int i = 0; i < n; i += (len << 1)) {
                for (int j = 0; j < len; j++) {
                    long u = a[i + j];
                    long v = a[i + j + len];
                    a[i + j] = u + v;
                    a[i + j + len] = u - v;
                }
            }
        }
        if (inverse) {
            // For exactness, n must divide all entries (true if you did forward then inverse).
            for (int i = 0; i < n; i++) a[i] /= n;
        }
    }

    // XOR convolution: c[k] = sum_{i xor j = k} a[i]*b[j]
    public static long[] xorConvolution(long[] a, long[] b) {
        int n = 1;
        int need = Math.max(a.length, b.length);
        while (n < need) n <<= 1;

        long[] fa = Arrays.copyOf(a, n);
        long[] fb = Arrays.copyOf(b, n);

        fwhtXor(fa, false);
        fwhtXor(fb, false);
        for (int i = 0; i < n; i++) fa[i] *= fb[i];
        fwhtXor(fa, true);
        return fa;
    }
}

