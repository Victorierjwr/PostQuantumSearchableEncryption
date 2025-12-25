package qut.edu.au.pqcsebaseopt.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleD.sampleD;

public class SampleDTest {

    @Test
    public void testSampleD() {
        int n = 3;
        BigInteger q = BigInteger.valueOf(17);

        // Define a simple lattice basis (3D)
        BigIntMatrix B = new BigIntMatrix(3, 3, q);
        // Row vectors:
        B.set(0, 0, 2); B.set(0, 1, 0); B.set(0, 2, 0);
        B.set(1, 0, 1); B.set(1, 1, 2); B.set(1, 2, 0);
        B.set(2, 0, 0); B.set(2, 1, 1); B.set(2, 2, 2);

        double s = 5.0;
        double[] c = {1.5, 1.5, 1.5}; // Target center
        int nParam = 17; // Security parameter

        System.out.println("Running SampleD...");

        // The return type is BigInteger[], safe from overflow
        BigInteger[] result = sampleD(B, s, c, nParam);

        System.out.println("Sampled Vector v: " + Arrays.toString(result));
    }
}
