import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class PredictProteinCoding {
    private static int k = 5;
    private static String genome;
    private static Set<Integer> genBankProteinLocation;
    private static Map[] scoringMaps;
    private static List[] frameSeparatedORFs;
    private static List<ORF> allORFs;

    public static void main(String[]args) throws FileNotFoundException {
        loadGenome();
        loadGenBank();
        identifyORFs();
        mergeAllORFs();

        List<ORF> shortORFs = getShortORFs(allORFs, k+1, 50);
        List<ORF> longORFs = getLongORFs(allORFs, 1400);

        Map<String, Integer>[] freqMaps = generateFrequencies(longORFs);
        scoringMaps = generateScoringMaps(longORFs, freqMaps[0], freqMaps[1]);

        for (ORF orf: allORFs) {
            orf.calcuateMMscore();
        }

        // PART 1
        partA();
        partBCD(shortORFs.size(), longORFs.size(), genBankProteinLocation.size());
        partE(freqMaps[0], freqMaps[1]);
        partF();

        // PART 2, 3, & 4
        //generatePlotDataLength();
        //generatePlotDataMMScore();

        // PART 5!!!
        // generateScatterplotData(shortORFs);
        //generateScatterplotData(longORFs);
        double[] shortMedians = getMedianValues(shortORFs);
        System.out.println();
        double[] longMedians = getMedianValues(longORFs);
        System.out.println();
        double[] statsAt20 = TPRandFPR(shortMedians, longMedians, .2);
        System.out.println("\nTPR with Line at 20%: " + statsAt20[0]);
        System.out.println("FPR with Line at 20%: " + statsAt20[1] + "\n");
        //System.out.println("\nSCATTERPLOT DATA (ALL):");
        //generateScatterplotData();


        //generateProfessorFlashbulb(shortORFs, longORFs);
        //generatePlotDataBoth(shortMedians, longMedians);
    }

    // length on x, MM scores on y

    // merges orfs from all three reading frames in the genome.
    public static void mergeAllORFs() {
        allORFs = new ArrayList(frameSeparatedORFs[0]);
        allORFs.addAll(frameSeparatedORFs[1]);
        allORFs.addAll(frameSeparatedORFs[2]);
        Collections.sort(allORFs);
    }

    // returns an array containing the median length and median MM Score for a List of ORFs
    public static double[] getMedianValues(List<ORF> orfs) {
        // rval [median of ORF lengths, median of ORFs score]
        // rval = [short length median, long length median, short MMscore median, long MMscore median];
        double[] rval = new double[2];
        Collections.sort(orfs, new Comparator<ORF>() {
            @Override
            public int compare(ORF o1, ORF o2) {
                return o1.length - o2.length;
            }
        });
        rval[0] = orfs.get(orfs.size()/2).length;

        Collections.sort(orfs, new Comparator<ORF>() {
            @Override
            public int compare(ORF o1, ORF o2) {
                return (int)(o1.MMscore) - (int)o2.MMscore;
            }
        });
        rval[1] = orfs.get(orfs.size()/2).MMscore;
        System.out.print(Arrays.toString(rval));
        return rval;
    }

    // calcuates TPR and FPR values based on a linear classifer created with the training medians
    public static double[] TPRandFPR(double[] shortMedians, double[] longMedians, double percentage) {
        double m = -1 / ((longMedians[1] - shortMedians[1]) / (longMedians[0] - shortMedians[0]));
        double x1 = shortMedians[0] + ((longMedians[0] - shortMedians[0]) * percentage);
        double y1 = shortMedians[1] + ((longMedians[1] - shortMedians[1]) * percentage);
        int TP = 0; int FP = 0; int FN = 0; int TN = 0;
        for(ORF orf: allORFs) {
            if (orf.MMscore > m * (orf.length - x1) + y1) {
                if (orf.genBankMatch) {
                    TP++;
                } else {
                    FP++;
                }
            } else {
                if (orf.genBankMatch) {
                    FN++;
                } else {
                    TN++;
                }
            }
        }
        return new double[]{(TP * 1.0) / (TP + FN), (FP * 1.0) / (FP + TN)};
    }

    // generates and x and y array representing plots to put on the ROC curve for the joint length and MMscore classifer
    public static void generatePlotDataBoth(double[] shortMedians, double[] longMedians) {
        double thresholdTPR = 0.0;
        double thresholdFPR = 0.0;
        double percentage = 0.0;
        double[] FPR = new double[200];
        double[] TPR = new double[200];
        for (int i = 0; i < 200; i++) {
            double[] results = TPRandFPR(shortMedians, longMedians, 0.005* i);
            TPR[i] = results[0];
            FPR[i] = results[1];

            if (TPR[i] > .8) {
                System.out.println(TPR[i]);
                percentage = 0.005 * i;
                thresholdTPR = TPR[i];
                thresholdFPR = FPR[i];
            }
        }
        System.out.println("ROC DATA WITH RESPECT TO BOTH LENGTH AND MM SCORE:");
        System.out.println("TPR = " + Arrays.toString(FPR) + ";");
        System.out.println("FPR = " + Arrays.toString(TPR) + ";");
        System.out.println("80% TPR THRESHOLD: percentage=" + percentage + ", TPR=" + thresholdTPR + ", FPR="+ thresholdFPR + "\n");
    }

    // generates and x and y array representing plots to put on the ROC curve for the length classifer
    public static void generatePlotDataLength() {
        int maxLength = -1;
        for(ORF orf: allORFs) {
            if (orf.length > maxLength) {
                maxLength = orf.length;
            }
        }
        double thresholdTPR = 0.0;
        double thresholdFPR = 0.0;
        int thresholdLength = 0;
        double[] FPR = new double[(maxLength - 6)/ 3];
        double[] TPR = new double[(maxLength - 6)/ 3];
        for (int i = 0; i < (maxLength - 6) / 3; i++) {
            int TP = 0; int FP = 0; int TN = 0; int FN = 0;
            for(ORF orf: allORFs) {
                if (orf.length > (i * 3) + 6) { // if length requirement is met
                    if (orf.genBankMatch) {
                        TP++;
                    } else {
                        FP++;
                    }
                } else {
                    if (orf.genBankMatch) {
                        FN++;
                    } else {
                        TN++;
                    }
                }
            }
            TPR[i] = (TP * 1.0) / (TP + FN);
            FPR[i] = (FP * 1.0) / (FP + TN);

            if (TPR[i] > .8) {
                thresholdLength = (i*3) + 6;
                thresholdTPR = TPR[i];
                thresholdFPR = FPR[i];
            }
        }
        System.out.println("ROC DATA WITH RESPECT TO LENGTH:");
        System.out.println("x = " + Arrays.toString(FPR) + ";");
        System.out.println("y = " + Arrays.toString(TPR) + ";");
        System.out.println("80% TPR THRESHOLD: Length=" + thresholdLength + ", TPR=" + thresholdTPR + ", FPR="+ thresholdFPR + "\n");

    }

    // generates and x and y array representing plots to put on the ROC curve for the MMscore classifer
    public static void generatePlotDataMMScore() {
        double MMmax = 0;
        double MMmin = 0;
        for(ORF orf: allORFs) {
            if (orf.MMscore > MMmax) {
                MMmax = orf.MMscore;
            }
            if (orf.MMscore < MMmin) {
                MMmin = orf.MMscore;
            }
        }
        double thresholdTPR = 0.0;
        double thresholdFPR = 0.0;
        double thresholdMMscore = 0.0;
        double[] FPR = new double[(int)(MMmax - MMmin)];
        double[] TPR = new double[(int)(MMmax - MMmin)];
        for (int i = 0; i < (int)(MMmax - MMmin); i++) {
            int TP = 0; int FP = 0; int TN = 0; int FN = 0;
            for(ORF orf: allORFs) {
                if (orf.MMscore > MMmin + i) { // if mmScore requirement is met
                    if (orf.genBankMatch) {
                        TP++;
                    } else {
                        FP++;
                    }
                } else {
                    if (orf.genBankMatch) {
                        FN++;
                    } else {
                        TN++;
                    }
                }
            }
            TPR[i] = (TP * 1.0) / (TP + FN);
            FPR[i] = (FP * 1.0) / (FP + TN);

            if (TPR[i] > .8) {
                thresholdMMscore = MMmin + i;
                thresholdTPR = TPR[i];
                thresholdFPR = FPR[i];
            }
        }
        System.out.println("ROC DATA WITH RESPECT TO MM SCORE:");
        System.out.println("x = " + Arrays.toString(FPR) + ";");

        System.out.println("y = " + Arrays.toString(TPR) + ";");
        System.out.println("80% TPR THRESHOLD: MMscore=" + thresholdMMscore + ", TPR=" + thresholdTPR + ", FPR="+ thresholdFPR + "\n");
    }

    // generates x and y array for all ORF data points
    public static void generateScatterplotData() {
        generateScatterplotData(allORFs);
    }

    // generates x and y array for List of ORFs passed in as a parameter
    public static void generateScatterplotData(List<ORF> orfs) {
        int proteins = 0;
        int nonproteins = 0;
        for(ORF orf: orfs) {
            if (orf.genBankMatch) {
                proteins++;
            } else {
                nonproteins++;
            }
        }
        int[] proteinLength = new int[proteins];
        int[] nonproteinLength = new int[nonproteins];
        double[] proteinMMscore = new double[proteins];
        double[] nonproteinMMscore = new double[nonproteins];

        int i1 = 0;
        int i2 = 0;
        for(int i = 0; i < orfs.size(); i++) {
            ORF temp = orfs.get(i);
            if (temp.genBankMatch) {
                proteinLength[i1] = temp.length;
                proteinMMscore[i1] = temp.MMscore;
                i1++;
            } else {
                nonproteinLength[i2] = temp.length;
                nonproteinMMscore[i2] = temp.MMscore;
                i2++;
            }
        }

        System.out.println("x1 = " + Arrays.toString(proteinLength) + ";");
        System.out.printf("y1 = [ %.2f", proteinMMscore[0]);
        for (int i = 1; i < proteinMMscore.length; i++) {
            System.out.printf(", %.2f", proteinMMscore[i]);
        }
        System.out.println("];");

        System.out.println("x2 = " + Arrays.toString(nonproteinLength) + ";");
        System.out.printf("y2 = [ %.2f", nonproteinMMscore[0]);
        for (int i = 1; i < nonproteinMMscore.length; i++) {
            System.out.printf(", %.2f", nonproteinMMscore[i]);
        }
        System.out.println("];");

    }

    // loads the genome file
    public static void loadGenome() throws FileNotFoundException {
        genome = "";
        int seq = 0;
        try (Scanner scanner = new Scanner(new File("GCF_000091665.1_ASM9166v1_genomic.txt"))) {
            while (scanner.hasNextLine()) {
                String ln = scanner.nextLine();
                if (ln.charAt(0) == '>') {
                    seq++;
                } else if (ln.charAt(0) != '>' && seq == 1) {
                    ln = ln.toUpperCase();
                    ln = ln.replaceAll("[^ACGT]", "T");
                    genome += ln;
                }
            }
        }
    }

    // loads the gen bank file and stores a map of end locations for positive scoring ORFs
    public static void loadGenBank() throws FileNotFoundException {
        genBankProteinLocation = new HashSet<>();
        int seq = 0;
        try (Scanner scanner = new Scanner(new File("GCF_000091665.1_ASM9166v1_genomic.gff"))) {
            while (scanner.hasNextLine()) {
                String ln = scanner.nextLine();
                if (ln.charAt(0) != '#') {
                    String[] fields = ln.split("\t"); // split whitespace && fields[1].equals("Protein Homology")
                    if (fields[0].equals("NC_000909.1") && fields[2].equals("CDS") && fields[6].equals("+")) {
                        genBankProteinLocation.add(Integer.parseInt(fields[4]));
                    }
                }
            }
        }
    }

    // scans genome and compiles list of ORFs
    public static void identifyORFs() {
        frameSeparatedORFs = new List[3];
        allORFs = new ArrayList<>();
        for (int frame = 0; frame <= 2; frame++) {
            List<ORF> tempORFs = new ArrayList<ORF>();
            int ORFstart = frame;
            //int ORFstop = -1;
            for (int i = frame; i + 3 < genome.length(); i += 3) {
                String nucleotide = genome.substring(i, i + 3);
                if (nucleotide.equals("TAA") || nucleotide.equals("TAG") || nucleotide.equals("TGA")) { // if stop codon
                    ORF temp = new ORF(ORFstart + 1, i);
                    tempORFs.add(temp);
                    allORFs.add(temp);
                    ORFstart = i + 3;
                } else if (i + 6 >= genome.length()) { // if this is the LAST nucleotide and its not a stop codon
                    ORF temp = new ORF(ORFstart + 1, i + 3);
                    tempORFs.add(temp);
                    allORFs.add(temp);
                    ORFstart = i + 3;
                }
            }
            frameSeparatedORFs[frame] = tempORFs;
        }
    }

    // prints a list of ORFs
    public static void printORFs(List<ORF> ORFs) {
        for(ORF orf: ORFs) {
            System.out.println(orf);
        }
    }

    // generates a frequency map for a training set of ORFs
    public static Map<String, Integer>[] generateFrequencies(List<ORF> trustedORFs) {
        Map<String, Integer> foreground = new HashMap<>();
        Map<String, Integer> background = new HashMap<>();
        Map<String, Integer>[] frequencyMaps = new Map[]{foreground, background};

        for (ORF orf: trustedORFs) {
            String seq = orf.sequence;
            for (int i = 0; i + k + 1 < seq.length(); i++) {
                if (foreground.get(seq.substring(i, i + k + 1)) == null) {
                    foreground.put(seq.substring(i, i + k + 1), 1);
                } else {
                    int count = foreground.get(seq.substring(i, i + k + 1));
                    foreground.put(seq.substring(i, i + k + 1), count += 1);
                }
            }
        }

        for (ORF orf: trustedORFs) {
            String rc = orf.revcomp;
            for (int i = 0; i + k + 1 < rc.length(); i++) {
                if (background.get(rc.substring(i, i + k + 1)) == null) {
                    background.put(rc.substring(i, i + k + 1), 1);
                } else {
                    int count = background.get(rc.substring(i, i + k + 1));
                    background.put(rc.substring(i, i + k + 1), count += 1);
                }
            }
        }
        return frequencyMaps;
    }

    // generates a scoring map based on frequencies
    public static Map<String, Double>[] generateScoringMaps(List<ORF> trustedORFs, Map<String, Integer> p, Map<String, Integer> q) {
        // log probability of each possible length k sequence occurring
        Map<String, Double> independentP = new HashMap<>();
        Map<String, Double> independentQ = new HashMap<>();

        // log probability of String[k+1] occuring immediately after String.substring(0, k + 1) occurs
        Map<String, Double> conditionalP = new HashMap<>();
        Map<String, Double> conditionalQ = new HashMap<>();

        String[] bases = new String[]{"A", "C", "G", "T"};
        String[] kLengthPermuations = generatePermutations(k);

        // compute unconditional log probabilities
        int denominator = (int) Math.pow(k, 4);
        for (ORF orf: trustedORFs) {
            denominator += orf.length - k + 1;
        }
        for (String permutation : kLengthPermuations) {
            for (int i = 0; i < 4; i++) {
                p.putIfAbsent(permutation + bases[i], 0);
                q.putIfAbsent(permutation + bases[i], 0);
            }
        }

        for (String permutation : kLengthPermuations) {
            int numeratorP = p.get(permutation + "A") + p.get(permutation + "C") + p.get(permutation + "G") + p.get(permutation + "T");
            int numeratorQ = q.get(permutation + "A") + q.get(permutation + "C") + q.get(permutation + "G") + q.get(permutation + "T");
            //double logprobP = Math.log(independentP.get(permutation) / total);
            //double logprobQ = Math.log(independentQ.get(permutation) / total);
            independentP.put(permutation, Math.log((numeratorP + 1) * 1.0 / (denominator + 1)));
            independentQ.put(permutation, Math.log((numeratorQ + 1) * 1.0 / (denominator + 1)));
        }

        // compute dependent log probabilities:
        for (String permutation : kLengthPermuations) {
            // foreground
            for (int i = 0; i < 4; i++) {
                p.putIfAbsent(permutation + bases[i], 0);
                q.putIfAbsent(permutation + bases[i], 0);
            }
            int countAp = 1 + p.get(permutation + "A");
            int countCp = 1 + p.get(permutation + "C");
            int countGp = 1 + p.get(permutation + "G");
            int countTp = 1 + p.get(permutation + "T");
            conditionalP.put(permutation + "A", Math.log(((countAp * 1.0) / (countAp + countCp + countGp + countTp))));
            conditionalP.put(permutation + "C", Math.log(((countCp * 1.0) / (countAp + countCp + countGp + countTp))));
            conditionalP.put(permutation + "G", Math.log(((countGp * 1.0) / (countAp + countCp + countGp + countTp))));
            conditionalP.put(permutation + "T", Math.log(((countTp * 1.0) / (countAp + countCp + countGp + countTp))));

            // background
            int countAq = 1 + q.get(permutation + "A");
            int countCq = 1 + q.get(permutation + "C");
            int countGq = 1 + q.get(permutation + "G");
            int countTq = 1 + q.get(permutation + "T");
            conditionalQ.put(permutation + "A", Math.log(((countAq * 1.0) / (countAq + countCq + countGq + countTq))));
            conditionalQ.put(permutation + "C", Math.log(((countCq * 1.0) / (countAq + countCq + countGq + countTq))));
            conditionalQ.put(permutation + "G", Math.log(((countGq * 1.0) / (countAq + countCq + countGq + countTq))));
            conditionalQ.put(permutation + "T", Math.log(((countTq * 1.0) / (countAq + countCq + countGq + countTq))));
        }

    return new Map[]{independentP, independentQ, conditionalP, conditionalQ};
    }

    // generates all the possible k length genetic permunations
    public static String[] generatePermutations(int length) {
        if (length == 1) {
            return new String[]{"A", "C", "G", "T"};
        } else {
            String[] builders = generatePermutations(length - 1);
            String[] rval = new String[(int) Math.pow(4, length)];
            for (int i = 0; i < builders.length; i++) {
                rval[(i * 4) + 0] = builders[i] + "A";
                rval[(i * 4) + 1] = builders[i] + "C";
                rval[(i * 4) + 2] = builders[i] + "G";
                rval[(i * 4) + 3] = builders[i] + "T";
            }
            return rval;
        }
    }

    // scans input list of ORFs for those that meet min length requirement
    public static List<ORF> getLongORFs(List<ORF> ORFs, int minLength) {
        List<ORF> longORFs = new ArrayList<>();
        for (ORF orf: ORFs) {
            if (orf.length > minLength) {
                longORFs.add(orf);
            }
        }
        return longORFs;
    }

    // scans input list of ORFs for those that meet min and max length requirement (min is so that we can exclude ORFs
    // to small to calcualte MM score for
    public static List<ORF> getShortORFs(List<ORF> ORFs, int minLength, int maxLength) {
        List<ORF> shortORFs = new ArrayList<>();
        for (ORF orf: ORFs) {
            if (orf.length >= minLength && orf.length < maxLength) {
                shortORFs.add(orf);
            }
        }
        return shortORFs;
    }

    // for each of the three reading frames, print a summary of the first and last reading frame in each
    // note: open reading frames  from each frame are stored as values in the ORFs variable.
    public static void partA() {
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("--------------------------     TASK 1A:     -----------------------------");
        System.out.println("-------------------------------------------------------------------------\n");

        for (int i = 0; i < 3; i++) {
            System.out.println("READING FRAME " + (i + 1) + ":");
            System.out.println("Total Number of ORFs:\t" + frameSeparatedORFs[i].size());
            System.out.println("Summary of First ORF:\t" + frameSeparatedORFs[i].get(0));
            System.out.println("Summary of Last ORF :\t" + (frameSeparatedORFs[i].get(frameSeparatedORFs[i].size() - 1)));
            System.out.println("\n");
        }

    }

    public static void partBCD(int shortOrfSize, int longOrfSize, int genBankSize) {
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("----------------------     TASK 1B, 1C & 1D:     ------------------------");
        System.out.println("-------------------------------------------------------------------------\n");

        System.out.println("SHORT ORFs (length < 50):\t" + shortOrfSize);
        System.out.println("LONG ORFs (length > 1400):\t" + longOrfSize);
        System.out.println("POS-STRAND CDSs in GENBANK:\t" + genBankSize);
        System.out.println("\n");
    }

    public static void partD() {
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("--------------------------     TASK 1D:     -----------------------------");
        System.out.println("-------------------------------------------------------------------------\n");

        System.out.println("PLUS STRAND CDSs in GENBANK:\t" + genBankProteinLocation.size());
        System.out.println("\n");
    }

    public static void partE(Map<String, Integer> p, Map<String, Integer> q) {
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("--------------------------     TASK 1E:     -----------------------------");
        System.out.println("-------------------------------------------------------------------------\n");

        String[] bases = new String[]{"A", "C", "G", "T"};

        System.out.println("P(T|AAGxy):\n");
        System.out.println(" \t" + String.format("%11s", ("A   ")) + String.format("%11s", ("C   ")) + String.format("%11s", ("G   ")) + String.format("%11s", ("T   ")));
        for (int i = 0; i < bases.length; i++) {
            System.out.print(bases[i] + "\t");
            for (int j = 0; j < bases.length; j++) {
                int Acounts = p.get("AAG" + bases[i] + bases[j] + "A");
                int Ccounts = p.get("AAG" + bases[i] + bases[j] + "C");
                int Gcounts = p.get("AAG" + bases[i] + bases[j] + "G");
                int Tcounts = p.get("AAG" + bases[i] + bases[j] + "T");
                System.out.printf("%11.9s", " " + ((Tcounts * 1.0) / (Acounts + Ccounts + Gcounts + Tcounts)));
                //System.out.println("P(T|AAG" + bases[i] + bases[j] + ") = " + );
            }
            System.out.println();
        }
        System.out.println("\n");
        System.out.println("Q(T|AAGxy):\n");
        System.out.println(" \t" + String.format("%11s", ("A   ")) + String.format("%11s", ("C   ")) + String.format("%11s", ("G   ")) + String.format("%11s", ("T   ")));
        for (int i = 0; i < bases.length; i++) {
            System.out.print(bases[i] + "\t");
            for (int j = 0; j < bases.length; j++) {
                int Acounts = q.get("AAG" + bases[i] + bases[j] + "A");
                int Ccounts = q.get("AAG" + bases[i] + bases[j] + "C");
                int Gcounts = q.get("AAG" + bases[i] + bases[j] + "G");
                int Tcounts = q.get("AAG" + bases[i] + bases[j] + "T");
                System.out.printf("%11.9s", " " + ((Tcounts * 1.0) / (Acounts + Ccounts + Gcounts + Tcounts)));
                //System.out.println("P(T|AAG" + bases[i] + bases[j] + ") = " + );
            }
            System.out.println();
        }
        System.out.println("\n");
    }

    public static void partF() {
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("--------------------------     TASK 1F:     -----------------------------");
        System.out.println("-------------------------------------------------------------------------\n");

        List<ORF> shortORFs = getShortORFs(allORFs, 6, 50);
        System.out.println("FIRST 5 SHORT ORFs SUMMARIES: ");
        for (int i = 0; i < 5; i++) {
            System.out.println(shortORFs.get(i));
        }
        System.out.println("\n");

        List<ORF> longORFs = getLongORFs(allORFs, 1400);
        System.out.println("FIRST 5 LONG ORFs SUMMARIES: ");
        for (int i = 0; i < 5; i++) {
            System.out.println(longORFs.get(i));
        }
        System.out.println("\n");
    }

    private static class ORF implements Comparable<ORF>{
        String sequence;
        String revcomp;
        int start;
        int stop;
        int length;
        boolean genBankMatch;
        double MMscore;

        boolean isRevComp;

        public ORF(int start, int stop) {
            this.sequence = genome.substring(start - 1, stop);
            this.revcomp = calculateReverseComplement(sequence);
            this.start = start;
            this.stop = stop;
            this.length = stop - start + 1;
            this.MMscore = 0.0;
            this.genBankMatch = genBankProteinLocation.contains(stop + 3);
        }

        // takes a sequence and returns reverse complement
        // e.g. AAGTTA...CAA => TTG...TAACTT
        public String calculateReverseComplement(String seq) {
            String rcomp = "";
            for (int i = 0; i < seq.length(); i++)
                rcomp = switch (seq.charAt(i)) {
                    case 'C' -> "G" + rcomp;
                    case 'G' -> "C" + rcomp;
                    case 'A' -> "T" + rcomp;
                    default -> "A" + rcomp;
                };
            return rcomp;
        }

        public void calcuateMMscore() {
            if (sequence.length() > k) {
                Map<String, Double> independentP = scoringMaps[0];
                Map<String, Double> independentQ = scoringMaps[1];
                Map<String, Double> conditionalP = scoringMaps[2];
                Map<String, Double> conditionalQ = scoringMaps[3];

                double pScore = independentP.get(this.sequence.substring(0, k));
                double qScore = independentQ.get(this.sequence.substring(0, k));
                for (int i = 1; i + k < sequence.length(); i++) {
                    pScore += conditionalP.get(this.sequence.substring(i, i + k + 1));
                    qScore += conditionalQ.get(this.sequence.substring(i, i + k + 1));
                }

                this.MMscore = pScore - qScore;
            }
        }

        public String toString() {
            return "start: " + String.format("%-10s", (start + ",")) + "stop: " + String.format("%-10s", (stop + ",")) + "length: " + String.format("%-9s", (length + ",")) +  "known match: " + String.format("%-9s", (genBankMatch + ",")) + "MM Score: " +String.format("%-9.6f", MMscore);
        }

        @Override
        public int compareTo(ORF o) {
            return start - o.start;
        }
    }
}
