import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

/**
 * Class: Sequences
 * 
 * Definition for list of sequences
 * 
 * @author gopi
 * 
 */
public class Sequences {
	// list of sequences
	ArrayList<Sequence> sequences = new ArrayList<Sequence>();

	// ORDER: A C G T
	int dna[] = new int[4];
	// ORDER: A R N D C E Q G H I L K M F P S T W Y V U O
	int aa[] = new int[20];
	String type = "";

	/*
	 * method to read the input sequence file
	 */
	void readSequences(String filename) {
		BufferedReader seqReader = null;
		try {
			seqReader = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			System.out.println("The sequences source file: " + filename
					+ ", was not found!");
		}
		String line = "";
		String seq = "";

		try {
			while ((line = seqReader.readLine()) != null) {
				if (line.startsWith(">")) {
					seq = "";
					continue;
				}

				if (line.equals("")) {
					sequences.add(new Sequence(seq.trim()));
					seq = "";
				}

				seq += line.trim();
				sequences.add(new Sequence(seq));
			}
		} catch (IOException e) {
			System.out.println("Could'nt read the sequence source file: "
					+ filename);
		}
		// get counts for each residue
		countSymbols();
	}

	/*
	 * method to get the counts for each residue of the sequences
	 */
	private void countSymbols() {
		// check validity of the sequences, DNA or AA
		if (!checkAll()) {
			System.out.println("All the sequences should be DNA or AA.");
			System.exit(0);
		}

		if (type.equals("dna")) {
			dna[0] = getCount("A");
			dna[1] = getCount("C");
			dna[2] = getCount("G");
			dna[3] = getCount("T");
		} else {
			aa[0] = getCount("A");
			aa[1] = getCount("R");
			aa[2] = getCount("N");
			aa[3] = getCount("D");
			aa[4] = getCount("C");
			aa[5] = getCount("E");
			aa[6] = getCount("Q");
			aa[7] = getCount("G");
			aa[8] = getCount("H");
			aa[9] = getCount("I");
			aa[10] = getCount("L");
			aa[11] = getCount("K");
			aa[12] = getCount("M");
			aa[13] = getCount("F");
			aa[14] = getCount("P");
			aa[15] = getCount("S");
			aa[16] = getCount("T");
			aa[17] = getCount("W");
			aa[18] = getCount("Y");
			aa[19] = getCount("V");
			// aa[20] = getCount("U");
			// aa[21] = getCount("O");
		}
	}

	/*
	 * method to obtain the counts of a given residue
	 */
	private int getCount(String symbol) {
		int count = 0;
		for (int index0 = 0; index0 < sequences.size(); index0++)
			for (int index = 0; index < sequences.get(index0).sequence.length(); index++) {
				if (sequences.get(index0).sequence.charAt(index) == symbol
						.charAt(0))
					count++;
			}
		return count;
	}

	/*
	 * method to check whether all the sequences are either DNA or AA
	 */
	private boolean checkAll() {
		for (int index = 1; index < sequences.size(); index++)
			if (!sequences.get(index - 1).type
					.equals(sequences.get(index).type)) {
				System.out.println("Check sequences: " + (index) + "-"
						+ (index + 1));
				return false;
			}
		type = sequences.get(0).type;
		return true;
	}

	/*
	 * method to print the sequences
	 */
	void printSequences() {
		System.out.println("Original Sequences:");
		for (int index = 0; index < sequences.size(); index++) {
			System.out.println(sequences.get(index));
		}
	}
}

/**
 * Class: Sequence
 * 
 * Definition of a single sequence
 * 
 * @author gopi
 * 
 */
class Sequence {
	// the sequence storage
	String sequence = "";
	// type, DNA or AA
	String type = "";
	// list of motifs
	ArrayList<Motif> motifSubStrings = new ArrayList<Motif>();

	// best index of the motif
	int bestIndex = 0;
	// highest log-2-odd index
	int topIndex = -1;
	double topLog2Odd = 0.0;

	// ORDER: A C G T
	int dna[] = new int[4];
	// ORDER: A R N D C E Q G H I L K M F P S T W Y V U O
	int aa[] = new int[20];

	// constructor
	public Sequence(String seqs) {
		sequence = seqs;
		// get the symbol counts to construct the residue count and frequency
		// matrices
		getSymbolCounts();
		// get the motifs of the sequence
		getMotifSubStrings();
		// get the initial random motif
		getRandomMotifIndex();
	}

	// method to select a random motif
	private void getRandomMotifIndex() {
		Random rand = new Random();
		// rand.nextInt(max - min + 1) + min
		bestIndex = rand.nextInt(motifSubStrings.size());
	}

	// method to extract motifs of the sequence
	private void getMotifSubStrings() {
		int count = sequence.length() - DiscoverMotifs.motifLength + 1;
		for (int index = 0; index < count; index++)
			motifSubStrings.add(new Motif(sequence.substring(index, index
					+ DiscoverMotifs.motifLength)));
	}

	// get the counts of each residue
	private void getSymbolCounts() {
		if (type(sequence).equals("dna")) {
			type = "dna";
			dna[0] = getCount(sequence, "A");
			dna[1] = getCount(sequence, "C");
			dna[2] = getCount(sequence, "G");
			dna[3] = getCount(sequence, "T");
		} else {
			type = "aa";
			aa[0] = getCount(sequence, "A");
			aa[1] = getCount(sequence, "R");
			aa[2] = getCount(sequence, "N");
			aa[3] = getCount(sequence, "D");
			aa[4] = getCount(sequence, "C");
			aa[5] = getCount(sequence, "E");
			aa[6] = getCount(sequence, "Q");
			aa[7] = getCount(sequence, "G");
			aa[8] = getCount(sequence, "H");
			aa[9] = getCount(sequence, "I");
			aa[10] = getCount(sequence, "L");
			aa[11] = getCount(sequence, "K");
			aa[12] = getCount(sequence, "M");
			aa[13] = getCount(sequence, "F");
			aa[14] = getCount(sequence, "P");
			aa[15] = getCount(sequence, "S");
			aa[16] = getCount(sequence, "T");
			aa[17] = getCount(sequence, "W");
			aa[18] = getCount(sequence, "Y");
			aa[19] = getCount(sequence, "V");
			// aa[20] = getCount(sequence, "U");
			// aa[21] = getCount(sequence, "O");
		}
	}

	// obtain the type of the sequence, DNA or AA
	private String type(String seqs) {
		for (int index = 0; index < seqs.length(); index++) {
			if (seqs.charAt(index) != 'A' && seqs.charAt(index) != 'C'
					&& seqs.charAt(index) != 'G' && seqs.charAt(index) != 'T')
				return "aa";
		}
		return "dna";
	}

	// get count of a given residue symbol
	private int getCount(String seq, String symbol) {
		int count = 0;
		for (int index = 0; index < seq.length(); index++) {
			if (seq.charAt(index) == symbol.charAt(0))
				count++;
		}
		return count;
	}

	// clear log-2-odds scores
	public void clearLogOdds() {
		for (int i = 0; i < motifSubStrings.size(); i++)
			motifSubStrings.get(i).log2odds = -1;
	}
}

/**
 * Class: Motif
 * 
 * Definition for the sequence motifs
 * @author gopi
 *
 */
class Motif {
	// the motif
	String subString = "";
	// the log-2-odds score of the motif
	double log2odds;

	public Motif(String str) {
		this.subString = str;
	}
}