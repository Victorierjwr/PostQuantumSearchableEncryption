package qut.edu.au.pqcsebase.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebase.algorithms.RandBasis;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;
import qut.edu.au.pqcsebase.algorithms.TrapGenFirst;

import java.math.BigInteger;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static qut.edu.au.pqcsebase.algorithms.RandBasis.randBasis;
import static qut.edu.au.pqcsebase.algorithms.RandBasis.randBasisPre;
import static qut.edu.au.pqcsebase.algorithms.SamplePre.computeGaussParam;
import static qut.edu.au.pqcsebase.algorithms.TrapGen.trapGen;

public class RandBasisTest {
    //@Test
    public long testRandBasis() {
        // === 1. Parameter Setup ===
        long n_long = 5;       // Lattice dimension n
        long q_long = 17;     // Modulus q (prime)
        BigInteger q = BigInteger.valueOf(q_long);

        // Gaussian parameter sigma
        // Note: RandBasis requires sigma to be sufficiently large (greater than the smoothing parameter of the basis).
        // If sigma is too small (e.g., 5.0), the sampled vectors will be too short and likely linearly dependent,
        double sigma; // Increased Gaussian parameter for better sampling
        long nParam = 5; // Security parameter for SampleZ

        System.out.println("=== 1. Generating valid Basis (A, S) using TrapGen ===");

        BigIntMatrix[] matrices = trapGen(n_long,q_long,1);
        BigIntMatrix A = matrices[0]; // Public matrix A
        BigIntMatrix S = matrices[1]; // Basis S such that A*S = 0 mod q
        sigma = computeGaussParam(S) * 20;
        // ---------------------------------------------------------------------

        System.out.println("Input S dimensions: " + S.getRowDimension() + " x " + S.getColumnDimension());

        // Verify input validity: A * S must be 0 mod q
        assertTrue("TrapGen Input Invalid: A*S must be 0", checkZero(A.multiplyQ(S)));
        System.out.println("=== 2. Running RandBasis to randomize basis S ===");
        long startTime = System.currentTimeMillis();
        BigIntMatrix S_randomized = randBasis(S, sigma, nParam);
        long endTime = System.currentTimeMillis();
        System.out.println("RandBasis completed in " + (endTime - startTime)/60000 + " minutes");
        System.out.println("Output S_randomized dimensions: " + S_randomized.getRowDimension() + " x " + S_randomized.getColumnDimension());

        System.out.println("=== 3. Verifying A * S_randomized = 0 mod q ===");
        assertTrue("RandBasis Output Invalid: A*S_randomized must be 0", checkZero(A.multiplyQ(S_randomized)));
        System.out.println("Verification successful: A * S_randomized = 0 mod q");

        return endTime - startTime;
    }

    //@Test
    public long testRandBasisPre(){
        long n_long = 5;       // Lattice dimension n
        long q_long = 17;     // Modulus q (prime)
        BigInteger q = BigInteger.valueOf(q_long);

        double sigma; // Increased Gaussian parameter for better sampling
        long nParam = 5; // Security parameter for SampleZ

        System.out.println("=== 1. Generating valid Basis (A, S) using TrapGen ===");

        BigIntMatrix[] matrices = trapGen(n_long,q_long,1);
        BigIntMatrix A = matrices[0]; // Public matrix A
        BigIntMatrix S = matrices[1]; // Basis S such that A*S = 0 mod q
        sigma = computeGaussParam(S) * 20;
        // ---------------------------------------------------------------------

        System.out.println("Input S dimensions: " + S.getRowDimension() + " x " + S.getColumnDimension());

        // Verify input validity: A * S must be 0 mod q
        assertTrue("TrapGen Input Invalid: A*S must be 0", checkZero(A.multiplyQ(S)));
        System.out.println("=== 2. Running RandBasisPre to randomize basis S ===");
        long startTime = System.currentTimeMillis();
        BigIntMatrix S_randomized = randBasisPre(S, sigma, nParam);
        long endTime = System.currentTimeMillis();
        System.out.println("RandBasisPre completed in " + (endTime - startTime)/60000 + " minutes");
        System.out.println("Output S_randomized dimensions: " + S_randomized.getRowDimension() + " x " + S_randomized.getColumnDimension());

        System.out.println("=== 3. Verifying A * S_randomized = 0 mod q ===");
        assertTrue("RandBasisPre Output Invalid: A*S_randomized must be 0", checkZero(A.multiplyQ(S_randomized)));
        System.out.println("Verification successful: A * S_randomized = 0 mod q");

        return endTime - startTime;
    }

    private boolean checkZero(BigIntMatrix bigIntMatrix) {
        for (int i = 0; i < bigIntMatrix.getRowDimension(); i++) {
            for (int j = 0; j < bigIntMatrix.getColumnDimension(); j++) {
                if (!bigIntMatrix.get(i, j).equals(BigInteger.ZERO)) {
                    return false;
                }
            }
        }
        return true;
    }

    @Test
    public void compareTime(){
        //long time1 = testRandBasis();
        long time2 = testRandBasisPre();
        //System.out.println("RandBasis time: " + time1 + " ms");
        System.out.println("RandBasisPre time: " + time2 + " ms");
    }

}
