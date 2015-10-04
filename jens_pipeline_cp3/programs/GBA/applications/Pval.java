package applications;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Pval {
	/**
	 * 
	 * @param args[0]: name of the file where the fasta sequences
	 * are there.
	 * args[1]: Number of protein sequence to work on
	 * args[2]: Start of the protein.
	 * args[3]: End of the protein sequence.
	 */
	public static void main(String []args)
	{
		//In the main we have to handle the
		//inputs
		/**************************************/
		BufferedReader inFile = null;
		try {
			inFile = new BufferedReader(new FileReader(args[0]));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//The protein count
		int proCount =Integer.parseInt(args[1]);
		int startPos = Integer.parseInt(args[2]);
		int endPos = Integer.parseInt(args[3]);
		String totalSeq = new String();
		//Local count of protein
		int localProCount =0;
		boolean reached = false;
		String line = null;
		try {
			while(null != (line = inFile.readLine()))
			{
				if(false == reached)
				{
					if(true == line.startsWith(">"))
					{
						localProCount++;
						if(localProCount == proCount)
						{
							reached=true;
						}
					}
					else
					{
						continue;
					}
				}
				else
				{
					if(true==line.equals(""))
					{
						break;
					}
					else
					{
						totalSeq = totalSeq.concat(line); 
					}
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			inFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//Now totalSeq contains the sequence that we require
		//Now we need to get the repeat
		String repeat = totalSeq.substring(startPos-1, endPos);
		//System.out.println("Sequence: " + totalSeq);
		//System.out.println("Repeat: " + repeat);
		//Now repeat contains the repeat sequence
		/**************************************/
		int customLen = repeat.length();
		
		Pval wa = new Pval();
		wa.getFrequency(totalSeq);
		
		//Initialize the complexity calculator
		ComplexityCalculator cc = new ComplexityCalculator();
		cc.initializeAlphabet();
		cc.initializeNewAlphabet();
	    cc.computeNewScoringMatrix();
	    cc.normalizeNewScoringMatrixPow();
		//Now create the vector of the repeats
		double total =0.0;
		double sqTotal =0.0;
		int count =0;
		for(int i =0; i < 500; i++)
		{
			String seq = wa.getRandomSeq(customLen);
			double value = cc.calculateNor2LetterEntropyWScoMatrix(seq);
			total += value;
			double sq = Math.pow(value,2.0);
			sqTotal += sq;
			count++;
		}
		double mean = total/ count;
		double variance = sqTotal/(count-1) - (Math.pow(mean, 2.0)*count)/(count-1);
		double lambda = Math.pow(mean, 3.0)/variance;
		
		//Now we have mean an lambda that we can use for calculating the
		//p value of the repeat.
		
		double complexity = cc.calculateNor2LetterEntropyWScoMatrix(repeat);
		double pvalue = wa.getPVale(complexity, mean, lambda);
		//System.out.println(pvalue);
		System.out.printf("%.3e\n", pvalue);	
			
	}
	
	
	public Pval()
	{
		randNum = new Random();
	}
	
	public void initCountFreq()
	{
		countFreq = new LinkedHashMap<Character, Integer>();
		
	    countFreq.put( 'A',0 );
	    countFreq.put( 'R',0 );
	    countFreq.put( 'N',0);
	    countFreq.put( 'D',0 );
	    countFreq.put( 'C',0 );
	    countFreq.put( 'Q',0 );
	    countFreq.put( 'E',0 );
	    countFreq.put( 'G',0 );
	    countFreq.put( 'H',0 ); 
	    countFreq.put( 'I',0 ); 
	    countFreq.put( 'L',0 );  
	    countFreq.put( 'K',0 );  
	    countFreq.put( 'M',0 );  
	    countFreq.put( 'F',0 ); 
	    countFreq.put( 'P',0 );  
	    countFreq.put( 'S',0 );
	    countFreq.put( 'T',0 );  
	    countFreq.put( 'W',0 );  
	    countFreq.put( 'Y',0 );  
	    countFreq.put( 'V',0 );
	}
	
	
	private double getPVale(double x, double mu, double lambda) {
		Runtime rt = Runtime.getRuntime();
		Process  pr = null;
		try {
			//System.out.println(x + " " + mu + " " + lambda);
			pr = rt.exec("./ProbIGD " + x + " " + mu + " " + lambda);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		BufferedReader runIn = new BufferedReader(new InputStreamReader(pr.getInputStream()));
		String pValStr = null;
		Double start=0.0, middle=0.0, end=0.0;
		try {
			pValStr = runIn.readLine();
			//Scan the line
			String regex = "^(\\S+)\\s+(\\S+)\\s+(\\S+)";
			Pattern pat = Pattern.compile(regex);
			Matcher mat = pat.matcher(pValStr);
			if(true == mat.find())
			{
				String first = mat.group(1);
				String second = mat.group(2);
				String third = mat.group(3);
				start = Double.parseDouble(first);
				middle = Double.parseDouble(second);
				end = Double.parseDouble(third);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		double pvalue = Math.exp(start) + Math.exp(middle + end);
		
		//double pValue = Double.parseDouble(pValStr);
		//System.out.println(pValue);
		try {
			runIn.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//return pValue;
		return pvalue;
	}

	
	private String getRandomSeq(int customLen) {
		String seqStr = new String();
		for(int i=0; i < customLen; i++)
		{
			double num = randNum.nextDouble();
			//System.out.println("Random Number: " + num);
			char c = getCharFromFreq(num);
			String tempStr = String.valueOf(c);
			seqStr = seqStr.concat(tempStr);
		}
		return seqStr;
	}

	private char getCharFromFreq(double num) {
		Set<Character> keys = charFreq.keySet();
		Iterator<Character> it=keys.iterator();
		while(false != it.hasNext())
		{
			char c = it.next();
			double tempNum = charFreq.get(c);
			double diff = tempNum-num;
			//System.out.println("Diff: " + diff +" Char: " + c);
			if(diff >=0)
			{
				//System.out.println("Returning: " + c);
				return c;
			}
		}
		
		return 0;
	}

	private void getFrequency(String totalFile) {
		
		//Now totalFile will contain the toatalFile
		//Calculate the frequencies of letter in it
		initCountFreq();
		char []tempChar = totalFile.toCharArray();
		for(int i=0; i < tempChar.length; i++)
		{
			int tempCount = countFreq.get(tempChar[i]);
			tempCount++;
			countFreq.put(tempChar[i], tempCount);
		}
		
		Set<Character> keys = countFreq.keySet();
		Iterator<Character> it = keys.iterator();
		charFreq = new LinkedHashMap<Character, Double>();
		double up=0;
		while(false != it.hasNext())
		{
			Character temp = it.next();
			if(false == it.hasNext())
			{
				//This is the last Character
				up =1;
				charFreq.put(temp, up);
				//System.out.println("Char: " + temp + " Freq: " + up);
			}
			else
			{
				double freq = (double)countFreq.get(temp)/(double)tempChar.length;
				up += freq;
				charFreq.put(temp, up);
				//System.out.println("Char: " + temp + " Freq: " + up);
			}
		}
	}

	Random randNum = null;
	String comFile = null;
	LinkedHashMap<Character, Integer> countFreq = null;
	LinkedHashMap<Character, Double> charFreq = null;
	

}
