package qut.edu.au.pqcsebaseopt.algorithmTest;

import org.junit.Test;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleZ.sampleZ;

public class SampleZTest {

    @Test
    public void testSampleZ() {
        double s = 10.0;
        double c = 5.5; // 中心在 5.5
        int n = 17;

        System.out.println("Sampling from D_Z, " + s + ", " + c);
        for(int i=0; i<10; i++) {
            System.out.print(sampleZ(s, c, n) + " ");
        }
    }
}
