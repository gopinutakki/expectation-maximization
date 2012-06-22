import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Class: DiscoverMotifs
 * 
 * Definition for discovering the motifs in the provided sequences.
 * 
 * @author gopi
 * 
 */
public class DiscoverMotifs {

	// list of sequences
	static Sequences seqs = null;
	// motif length
	static int motifLength = 6;
	// initial random alignment count
	static int initialRandomCount = 50;
	// total iterations for each initial random alignment
	static int iterationsCount = 500;
	// input sequence file name
	static String inputFile = null;
	// pseudocount to prevent zero values
	static int pseudocount = 1;

	// matrices/tables to store counts and frequencies of residues
	static Residue counts = new Residue(1);
	static Residue freqs = new Residue(1);
	// output file writer
	static BufferedWriter bw = null;

	/*
	 * begins here... method to read the sequences, run configuration and start
	 * expectation.
	 */
	public static void main(String args[]) throws IOException {

		bw = new BufferedWriter(new FileWriter("output.txt"));
		System.out.println("Writing the results to 'output.txt'");
		try {
			readConfig();
		} catch (IOException e) {
			System.out
					.println("Could not read 'run.config' file for run configurations.");
		}
		seqs = new Sequences();
		seqs.readSequences(inputFile);

		for (int looper = 0; looper < initialRandomCount; looper++) {
			bw.write("\n\nInitial Random Alignment Run Number: " + (looper + 1)
					+ "\n");
			expectation();
			System.out.print(".");
		}

		printTopResults();
		bw.flush();
		bw.close();
		System.out.println("\nDONE!");
	}

	/*
	 * method to print the top log-2-odds scores
	 */
	private static void printTopResults() throws IOException {
		bw.write("\n\nTop Log-base-2 Odd scores (max log-2-odds of each sequence, from all runs together):\n");
		bw.write("Table columns in the order of: Sequence Number, Best Motif Location, Best Motif, Max-Log2Odd\n");
		bw.write("--------------------------------------------------------------------------------------------\n");
		for (int loop = 0; loop < seqs.sequences.size(); loop++) {
			bw.write((loop + 1)
					+ "\t\t"
					+ (seqs.sequences.get(loop).topIndex + 1)
					+ "\t\t"
					+ seqs.sequences.get(loop).motifSubStrings
							.get(seqs.sequences.get(loop).topIndex).subString
					+ "\t\t" + seqs.sequences.get(loop).topLog2Odd + "\n");
		}
	}

	/*
	 * method to initiate the expectation stage.
	 */
	private static void expectation() throws IOException {

		bw.write("Table columns in the order of: Sequence Number, Best Motif Location, Best Motif, Max-Log2Odd\n");
		bw.write("--------------------------------------------------------------------------------------------\n");
		// for all sequences, calculate residue counts and frequencies.
		for (int loop = 0; loop < seqs.sequences.size(); loop++) {
			for (int count = 0; count < iterationsCount; count++) {
				counts = new Residue(motifLength + 1);
				freqs = new Residue(motifLength + 1);
				// seqs.sequences.get(loop).clearLogOdds();
				if (seqs.sequences.get(0).type.equals("dna")) {
					expectationDNA(loop);
				} else
					expectationAA(loop);
			}
			// write the result to the output file.
			bw.write((loop + 1)
					+ "\t\t"
					+ (seqs.sequences.get(loop).bestIndex + 1)
					+ "\t\t"
					+ seqs.sequences.get(loop).motifSubStrings
							.get(seqs.sequences.get(loop).bestIndex).subString
					+ "\t\t"
					+ seqs.sequences.get(loop).motifSubStrings
							.get(seqs.sequences.get(loop).bestIndex).log2odds
					+ "\n");

			// printFreq();
		}

	}

	/*
	 * method to perform expectation stage for Amino Acids.
	 */
	private static void expectationAA(int loop) {
		String randomMotif = "";
		int location = -2;
		int totalPerSymbol = 0;
		int totalBackground = 0;

		// calculate counts for the residues, also background scores.
		for (int i = 0; i < 20; i++) {
			for (int j = 1; j < motifLength + 1; j++) {
				counts.aa[i][j] = 0;
				for (int index1 = 0; index1 < seqs.sequences.size(); index1++) {
					location = seqs.sequences.get(index1).bestIndex;
					// System.out.println(index1 + "--" + location);
					randomMotif = seqs.sequences.get(index1).motifSubStrings
							.get(location).subString;

					if (i == 0 && randomMotif.charAt(j - 1) == 'A')
						counts.aa[i][j]++;
					if (i == 1 && randomMotif.charAt(j - 1) == 'R')
						counts.aa[i][j]++;
					if (i == 2 && randomMotif.charAt(j - 1) == 'N')
						counts.aa[i][j]++;
					if (i == 3 && randomMotif.charAt(j - 1) == 'D')
						counts.aa[i][j]++;
					if (i == 4 && randomMotif.charAt(j - 1) == 'C')
						counts.aa[i][j]++;
					if (i == 5 && randomMotif.charAt(j - 1) == 'E')
						counts.aa[i][j]++;
					if (i == 6 && randomMotif.charAt(j - 1) == 'Q')
						counts.aa[i][j]++;
					if (i == 7 && randomMotif.charAt(j - 1) == 'G')
						counts.aa[i][j]++;
					if (i == 8 && randomMotif.charAt(j - 1) == 'H')
						counts.aa[i][j]++;
					if (i == 9 && randomMotif.charAt(j - 1) == 'I')
						counts.aa[i][j]++;
					if (i == 10 && randomMotif.charAt(j - 1) == 'L')
						counts.aa[i][j]++;
					if (i == 11 && randomMotif.charAt(j - 1) == 'K')
						counts.aa[i][j]++;
					if (i == 12 && randomMotif.charAt(j - 1) == 'M')
						counts.aa[i][j]++;
					if (i == 13 && randomMotif.charAt(j - 1) == 'F')
						counts.aa[i][j]++;
					if (i == 14 && randomMotif.charAt(j - 1) == 'P')
						counts.aa[i][j]++;
					if (i == 15 && randomMotif.charAt(j - 1) == 'S')
						counts.aa[i][j]++;
					if (i == 16 && randomMotif.charAt(j - 1) == 'T')
						counts.aa[i][j]++;
					if (i == 17 && randomMotif.charAt(j - 1) == 'W')
						counts.aa[i][j]++;
					if (i == 18 && randomMotif.charAt(j - 1) == 'Y')
						counts.aa[i][j]++;
					if (i == 19 && randomMotif.charAt(j - 1) == 'V')
						counts.aa[i][j]++;
					// if (i == 20 && randomMotif.charAt(j - 1) == 'U')
					// counts.aa[i][j]++;
					// if (i == 21 && randomMotif.charAt(j - 1) == 'O')
					// counts.aa[i][j]++;
				}
			}

			// add the pseudocounts
			for (int j = 1; j < motifLength + 1; j++) {
				totalPerSymbol += counts.aa[i][j];
				counts.aa[i][j] += pseudocount;
			}
			// get the background scores
			counts.aa[i][0] = seqs.aa[i] - totalPerSymbol;
			totalBackground += counts.aa[i][0];
		}

		// calculate frequencies
		for (int i = 0; i < 20; i++) {
			for (int j = 1; j < motifLength + 1; j++) {
				freqs.aa[i][j] = round((counts.aa[i][j])
						/ (seqs.sequences.size() - 1 + (pseudocount * 20)));
			}
			freqs.aa[i][0] = round((counts.aa[i][0] + 1)
					/ (totalBackground + (pseudocount * 20)));
		}

		// start maximization
		seqs.sequences.get(loop).bestIndex = maximizationAA(loop);
		// printFreq();
	}

	/*
	 * method to perform expectation stage for DNA.
	 */
	private static void expectationDNA(int loop) {

		String randomMotif = "";
		int location = -1;
		int totalPerSymbol = 0;
		int totalBackground = 0;

		// calculate counts.
		for (int i = 0; i < 4; i++) {
			for (int j = 1; j < motifLength + 1; j++) {
				counts.dna[i][j] = 0;
				for (int index1 = 0; index1 < seqs.sequences.size(); index1++) {
					location = seqs.sequences.get(index1).bestIndex;
					randomMotif = seqs.sequences.get(index1).motifSubStrings
							.get(location).subString;

					if (i == 0 && randomMotif.charAt(j - 1) == 'A')
						counts.dna[i][j]++;
					if (i == 1 && randomMotif.charAt(j - 1) == 'C')
						counts.dna[i][j]++;
					if (i == 2 && randomMotif.charAt(j - 1) == 'G')
						counts.dna[i][j]++;
					if (i == 3 && randomMotif.charAt(j - 1) == 'T')
						counts.dna[i][j]++;
				}
			}
			for (int j = 1; j < motifLength + 1; j++) {
				totalPerSymbol += counts.dna[i][j];
				// apply pseudocount.
				counts.dna[i][j] += pseudocount;
			}
			counts.dna[i][0] = seqs.dna[i] - totalPerSymbol;
			totalBackground += counts.dna[i][0];
		}

		// calculate frequencies.
		for (int i = 0; i < 4; i++) {
			for (int j = 1; j < motifLength + 1; j++) {
				freqs.dna[i][j] = round((counts.dna[i][j])
						/ (seqs.sequences.size() - 1 + (pseudocount * 4)));
			}
			freqs.dna[i][0] = round((counts.dna[i][0] + 1)
					/ (totalBackground + (pseudocount * 4)));
		}

		// start maximization
		seqs.sequences.get(loop).bestIndex = maximizationDNA(loop);
		// printFreq();
	}

	/*
	 * method to perform maximization for DNA sequences
	 */
	private static int maximizationDNA(int loopIndex) {
		double nume = 1, denom = 1;
		// calculate max log-2-odd scores for all the motifs of a given
		// sequence.
		for (int loop = 0; loop < seqs.sequences.get(loopIndex).motifSubStrings
				.size(); loop++) {
			String tempMotif = seqs.sequences.get(loopIndex).motifSubStrings
					.get(loop).subString;
			for (int index = 0; index < motifLength; index++) {
				nume *= freqs.dna[getDNASymbolIndex(tempMotif.charAt(index))][index];
				denom *= freqs.dna[getDNASymbolIndex(tempMotif.charAt(index))][0];
			}
			nume = round(nume);
			denom = round(denom);
			seqs.sequences.get(loopIndex).motifSubStrings.get(loop).log2odds = round(Math
					.log(nume / denom) / Math.log(2.0));

			nume = denom = 1;
		}

		// find the max log-2-odd scores
		double tempVal = -1;
		int maxIndex = -1;
		for (int loop = 0; loop < seqs.sequences.get(loopIndex).motifSubStrings
				.size(); loop++) {
			if (tempVal < seqs.sequences.get(loopIndex).motifSubStrings
					.get(loop).log2odds) {
				tempVal = seqs.sequences.get(loopIndex).motifSubStrings
						.get(loop).log2odds;
				maxIndex = loop;
			}
		}

		// keep track of the highest log-2-odd scores
		if (tempVal > seqs.sequences.get(loopIndex).topLog2Odd) {
			seqs.sequences.get(loopIndex).topLog2Odd = tempVal;
			seqs.sequences.get(loopIndex).topIndex = maxIndex;
			// System.out.println(seqs.sequences.get(loopIndex).topLog2Odd);
		}
		return maxIndex;
	}

	/*
	 * method to perform maximization for AA sequences
	 */
	private static int maximizationAA(int loopIndex) {
		double nume = 1, denom = 1;
		// calculate max log-2-odd scores for all the motifs of a given
		// sequence.
		for (int loop = 0; loop < seqs.sequences.get(loopIndex).motifSubStrings
				.size(); loop++) {
			String tempMotif = seqs.sequences.get(loopIndex).motifSubStrings
					.get(loop).subString;
			for (int index = 0; index < motifLength; index++) {
				nume *= freqs.aa[getAASymbolIndex(tempMotif.charAt(index))][index];
				denom *= freqs.aa[getAASymbolIndex(tempMotif.charAt(index))][0];
			}
			nume = round(nume);
			denom = round(denom);
			seqs.sequences.get(loopIndex).motifSubStrings.get(loop).log2odds = round(Math
					.log(nume / denom) / Math.log(2.0));

			nume = denom = 1;
		}

		// find the max log-2-odd scores
		double tempVal = -1;
		int maxIndex = -1;
		for (int loop = 0; loop < seqs.sequences.get(loopIndex).motifSubStrings
				.size(); loop++) {
			if (tempVal < seqs.sequences.get(loopIndex).motifSubStrings
					.get(loop).log2odds) {
				tempVal = seqs.sequences.get(loopIndex).motifSubStrings
						.get(loop).log2odds;
				maxIndex = loop;
			}
		}

		// keep track of the highest log-2-odd scores
		if (tempVal > seqs.sequences.get(loopIndex).topLog2Odd) {
			seqs.sequences.get(loopIndex).topLog2Odd = tempVal;
			seqs.sequences.get(loopIndex).topIndex = maxIndex;
			// System.out.println(seqs.sequences.get(loopIndex).topLog2Odd);
		}
		return maxIndex;
	}

	/*
	 * method to help populate the residue counts, frequencies matrices
	 */
	private static int getDNASymbolIndex(char ch) {
		if (ch == 'A')
			return 0;
		if (ch == 'C')
			return 1;
		if (ch == 'G')
			return 2;
		return 3;
	}

	/*
	 * method to help populate the residue counts, frequencies matrices
	 */
	private static int getAASymbolIndex(char ch) {
		// A R N D C E Q G H I L K M F P S T W Y V U O
		if (ch == 'A')
			return 0;
		if (ch == 'R')
			return 1;
		if (ch == 'N')
			return 2;
		if (ch == 'D')
			return 3;
		if (ch == 'C')
			return 4;
		if (ch == 'E')
			return 5;
		if (ch == 'Q')
			return 6;
		if (ch == 'G')
			return 7;
		if (ch == 'H')
			return 8;
		if (ch == 'I')
			return 9;
		if (ch == 'L')
			return 10;
		if (ch == 'K')
			return 11;
		if (ch == 'M')
			return 12;
		if (ch == 'F')
			return 13;
		if (ch == 'P')
			return 14;
		if (ch == 'S')
			return 15;
		if (ch == 'T')
			return 16;
		if (ch == 'W')
			return 17;
		if (ch == 'Y')
			return 18;
		// if (ch == 'V')
		return 19;
		// if (ch == 'U')
		// return 20;
		// if (ch == 'O')
		// return 21;
	}

	/*
	 * method to read the run configuration file
	 */
	private static void readConfig() throws IOException {
		BufferedReader br = new BufferedReader(new FileReader("run.config"));
		String line = "";

		while ((line = br.readLine()) != null) {

			if (line.startsWith("#"))
				continue;
			if (line.startsWith("#Comments"))
				break;

			if (line.startsWith("motif-count")) {
				String[] opts1 = line.split(":");
				if (opts1.length == 2) {
					motifLength = Integer.parseInt(opts1[1].trim());
				}
			} else if (line.startsWith("initial-random-count")) {
				String[] opts1 = line.split(":");
				if (opts1.length == 2) {
					initialRandomCount = Integer.parseInt(opts1[1].trim());
				}
			} else if (line.startsWith("iterations-count")) {
				String[] opts1 = line.split(":");
				if (opts1.length == 2) {
					iterationsCount = Integer.parseInt(opts1[1].trim());
				}
			} else if (line.startsWith("input-seq-file")) {
				String[] opts1 = line.split(":");
				if (opts1.length == 2) {
					inputFile = opts1[1].trim();
				}
			}
		}
	}

	/*
	 * method for decimal rounding
	 */
	private static double round(double d) {
		// try {
		// DecimalFormat twoDForm = new DecimalFormat("#.######");
		// return Double.valueOf(twoDForm.format(d));
		// } catch (NumberFormatException nfe) {
		// return 0.0;
		// }
		return d;
	}
}
