package qut.edu.au.pqcsebaseopt.projectsTest.IBEKSTEST;

import org.junit.Test;
import qut.edu.au.pqcsebaseopt.projects.IBEKS.Setup;
import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.*;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;

import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Encrypt.encrypt;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.KeyGen.keyGen;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Test.test;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Trapdoor.trapdoor;

public class SetupTest {

    @Test
    public void testSetup() {
        Setup setup = new Setup();
        System.out.println("Setup instance created: " + setup);
    }

    @Test
    public void keyGenTest(){
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());
    }

    @Test
    public void encryptTest(){
        String keyword = "engineer";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);

    }

    @Test
    public void trapdoorTest(){
        String keyword = "engineer";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);

        //---------------------------Trapdoor-----------------------------------
        String[] keywords = new String[]{"doctor", "engineer", "teacher"};
        TW tw = trapdoor(setup,keywords,keyPair.sk_id,identity);
        System.out.println("Generated Trapdoor: " + tw);

    }

    @Test
    public void testAll(){

        long startTime = System.currentTimeMillis();
        long endTime;

        String keyword = "engine";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Setup-----------------------------------");
        long setupTime = endTime - startTime;
        System.out.println("Setup Time: " + (endTime - startTime) + " ms");

        startTime = System.currentTimeMillis();
        KeyPair keyPair = keyGen(setup,msk,identity);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------KeyGen-----------------------------------");
        long keyGenTime = endTime - startTime;
        System.out.println("KeyGen Time: " + (endTime - startTime) + " ms");

        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        startTime = System.currentTimeMillis();
        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Encrypt-----------------------------------");
        long encryptTime = endTime - startTime;
        System.out.println("Encrypt Time: " + (endTime - startTime) + " ms");

        //---------------------------Trapdoor-----------------------------------
        String[] keywords = new String[]{"doctor", "engineer", "teacher"};

        startTime = System.currentTimeMillis();
        TW tw = trapdoor(setup,keywords,keyPair.sk_id,identity);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Trapdoor-----------------------------------");
        long trapdoorTime = endTime - startTime;
        System.out.println("Trapdoor Time: " + (endTime - startTime) + " ms");
        System.out.println("Generated Trapdoor: " + tw);

        startTime = System.currentTimeMillis();
        Boolean flag = test(ct,tw);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Test-----------------------------------");
        long testTime = endTime - startTime;
        System.out.println("Test Time: " + (endTime - startTime) + " ms");

        System.out.println("Test result: " + flag);

        System.out.println("----------------------------------------------------------");
        System.out.println("Summary of Times:");
        System.out.println("Setup Time: " + setupTime + " ms");
        System.out.println("KeyGen Time: " + keyGenTime + " ms");
        System.out.println("Encrypt Time: " + encryptTime + " ms");
        System.out.println("Trapdoor Time: " + trapdoorTime + " ms");
        System.out.println("Test Time: " + testTime + " ms");
        System.out.println("----------------------------------------------------------");

    }

    @Test
    public void hashTest(){
        Setup setup = new Setup();
        String testString = "123456789";
        System.out.println("Hashing string to matrix:");
        BigIntMatrix hashMatrix = setup.hash1(testString);
        System.out.println("Input String: " + testString);
        System.out.println("Hashed Matrix: " + hashMatrix);

        String keyword = "engineer";
        System.out.println("Hashing string to value:");
        BigInteger hashValue = setup.hash2(keyword);
        System.out.println("Input String: " + keyword);
        System.out.println("Hashed Value: " + hashValue);
    }
}
