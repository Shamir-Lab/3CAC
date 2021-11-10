
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class Algorithm3CAC {
	
	public static int binarySearch(ArrayList<Integer> list, int a){// return the index if exists, and return -(insertPos+1) if doesn't exist
		int l=0, r=list.size() - 1;
        while (l <= r) {
            int m = (l + r)/2;
            if (list.get(m) == a){
            	return m;
            }
            else if (list.get(m) < a){
            	l = m + 1;
            }
            else{
            	r = m - 1;
            }
        }
        return -(l+1);
	}
	
	public static int binaryInsert(ArrayList<Integer> list, int a){
		int l=0, r=list.size() - 1;
        while (l <= r) {
            int m = (l + r)/2;
            if (list.get(m) == a){
            	return -1;
            }
            else if (list.get(m) < a){
            	l = m + 1;
            }
            else{
            	r = m - 1;
            }
        }
        return l;
	}
	
	public static int getContigID(String contig, boolean flyeAssembly){
		if(flyeAssembly){
			return Integer.parseInt(contig.substring(7));
		}
		else{
			String[] temp=contig.trim().split("_");
			return Integer.parseInt(temp[1]);
		}
	}
	
	public static int getEdgeID(String edge, boolean flyeAssembly){
		if(flyeAssembly){
			return Integer.parseInt(edge.substring(5));
		}
		else{
			return Integer.parseInt(edge);
		}
	}
	
	public static int flyeContigsToEdgesInGraph(String file_assembly_info, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges) {
		try{
			//read edges of each contig
			int max_edge_ID=0;
			Scanner in=new Scanner(new File(file_assembly_info)); in.nextLine();  
			String[] readLine, edgeInPath;
			int contig_ID=-1, insert_pos=-1, edge_ID=-1;
			while(in.hasNextLine()) {
				readLine=in.nextLine().trim().split("[\\p{Space}]+");
				edgeInPath=readLine[7].trim().split(",");
				ArrayList<Integer> edgesInCurContig=new ArrayList<Integer>();
				for(int i=0;i<edgeInPath.length;i++) {
					if(!edgeInPath[i].equalsIgnoreCase("*")) {
						edge_ID = Math.abs(Integer.parseInt(edgeInPath[i]));
						edgesInCurContig.add(edge_ID);
						if(max_edge_ID < edge_ID){
							max_edge_ID = edge_ID;
						}
					}
				} 
				
				contig_ID=Integer.parseInt(readLine[0].substring(7));
				insert_pos=binaryInsert(contigs, contig_ID);
				contigs.add(insert_pos, contig_ID);
				contigEdges.add(insert_pos, edgesInCurContig);
			}
			in.close();
			System.out.println(contigs.size()+"  contigs contained in "+file_assembly_info+" !");
			return max_edge_ID;
			
		}catch(Exception e){
			
		}
		return -1;
	}	
	
	public static int spadesContigsToEdgesInGraph(String file_assembly_info, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges) {
		try{
			//read edges of each contig
			int max_edge_ID=0;
			Scanner in=new Scanner(new File(file_assembly_info));   
			ArrayList<Integer> edgesInCurContig=new ArrayList<Integer>();
			int contig_ID=-1, edge_ID=-1, insert_pos=-1, length=-1; boolean readNode=false;
			String readLine; String[] nodeLabel, edgeInPath;
			
			while(in.hasNextLine()) {
				readLine=in.nextLine();
				if(readLine.length()>4 && readLine.substring(0, 4).equalsIgnoreCase("NODE")){//read contig ID
					length=readLine.length();
					if(readLine.charAt(length-1)=='\''){
						readNode=false;
					}
					else{
						readNode=true;
						if(edgesInCurContig.size()>0){
							ArrayList<Integer> temp=new ArrayList<Integer>();
							for(int i=0;i<edgesInCurContig.size();i++){
								temp.add(edgesInCurContig.get(i));
							}
							contigEdges.add(insert_pos, temp);
							edgesInCurContig.clear();
						}
						
						nodeLabel=readLine.trim().split("_");
						contig_ID=Integer.parseInt(nodeLabel[1]);
						insert_pos=binaryInsert(contigs, contig_ID);
						contigs.add(insert_pos, contig_ID);
						
					}
				}
				else{
					if(readNode){//read the path if this Node was add to contigs
						length=readLine.length();
						if(readLine.charAt(length-1)==';'){
							readLine=readLine.substring(0, length-1);
						}
						edgeInPath=readLine.trim().split(",");
						for(int i=0;i<edgeInPath.length;i++){
							length=edgeInPath[i].length();
							edge_ID=Integer.parseInt(edgeInPath[i].substring(0, length-1));
							edgesInCurContig.add(edge_ID);
							
							if(max_edge_ID < edge_ID){
								max_edge_ID = edge_ID;
							}
						}
					}
				}
			}
			in.close();
			
			if(edgesInCurContig.size()>0){
				contigEdges.add(insert_pos, edgesInCurContig);
			}
			System.out.println(contigs.size()+"  contigs contained in "+file_assembly_info+" !");
			return max_edge_ID;
			
		}catch(Exception e){
			
		}
		return -1;
	}	
		
	public static void contigCategoryToEdges(String file_initial, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges, boolean flyeAssembly,
			int[] edgesCate){		
		
		try{
			Scanner in=new Scanner(new File(file_initial)); 
			String[] readLine; int contig_ID=-1, contig_cate=-1, contig_pos=-1, edge_ID=-1, edge_insert_pos=-1;
			
			while(in.hasNextLine()) {
				readLine = in.nextLine().trim().split(",");
				contig_ID = getContigID(readLine[0], flyeAssembly);
				
				if(readLine[1].equalsIgnoreCase("phage")) {
					contig_cate=1;	
				}
				else if(readLine[1].equalsIgnoreCase("plasmid")) {
					contig_cate=2;
				}
				else if(readLine[1].equalsIgnoreCase("chromosome")) {
					contig_cate=3;
				}
				else {
					contig_cate=4;
				}
				
				//assign category of contigs to edges
				if(contig_cate !=4){
					contig_pos = binarySearch(contigs, contig_ID);
					for(int i=0;i<contigEdges.get(contig_pos).size();i++) {
						edge_ID=contigEdges.get(contig_pos).get(i)-1;
						
						if(edgesCate[edge_ID] == 4){
							edgesCate[edge_ID]=contig_cate;
						}
						else if(edgesCate[edge_ID] != contig_cate){
							edgesCate[edge_ID]=-1; //edge with conflict categories as it was contained in multiple contigs with different categories
						}
					}
				}
			}
			in.close(); 
			
			for(int i=0; i<edgesCate.length; i++){
				if(edgesCate[i]==-1){ //set conflict edges to be uncertain
					edgesCate[i]=4; 
				}
			}
			System.out.println("Edges contained in contigs are assigned with contig categories ! ");
			
		}catch(Exception e){
			
		}
	}
	
	public static void readAssemblyGraph(String file_assembly_graph, boolean flyeAssembly, int[] edgesLength, ArrayList<ArrayList<Integer>> edgesAdj) {
		try{
			System.out.println("Reading assembly graph "+file_assembly_graph+" !");
			
			Scanner in=new Scanner(new File(file_assembly_graph));   
			if(flyeAssembly){
				in.nextLine();
			}
			
			String[] readLine; boolean end=false; 
			int edge_ID=-1, insert_pos=-1, preEdge_ID=-1, sufEdge_ID=-1;
			while(in.hasNextLine() && !end) {
				readLine=in.nextLine().trim().split("[\\p{Space}]+");
				if(readLine[0].equalsIgnoreCase("S")) {
					edge_ID=getEdgeID(readLine[1], flyeAssembly)-1;
					edgesLength[edge_ID] = readLine[2].length();
				}
				else if(readLine[0].equalsIgnoreCase("L")) {
					if(!readLine[1].equalsIgnoreCase(readLine[3])){
						preEdge_ID = getEdgeID(readLine[1], flyeAssembly)-1;
						sufEdge_ID = getEdgeID(readLine[3], flyeAssembly)-1;
						
						insert_pos=binarySearch(edgesAdj.get(preEdge_ID), sufEdge_ID);
						if(insert_pos < 0){
							insert_pos = Math.abs(insert_pos)-1;
							edgesAdj.get(preEdge_ID).add(insert_pos, sufEdge_ID);
							
							insert_pos=binarySearch(edgesAdj.get(sufEdge_ID), preEdge_ID);
							insert_pos = Math.abs(insert_pos)-1;
							edgesAdj.get(sufEdge_ID).add(insert_pos, preEdge_ID);
						}
					}
				}
				else {
					end=true;
				}
			}
			in.close();
			System.out.println("End Reading assembly graph !");
			
		}catch(Exception e){
			
		}
	}

	public static void correction(ArrayList<ArrayList<Integer>> edgesAdj, int[] edgesCate){
		
		int curEdgeCate=-1, curEdgeAdj=-1, adjCate=-1, correctEdges=0; 
		int num_support=0, num_same=0, loop_num=0;  
		boolean multiAdjCates=false, correct=true; 
		
		while(correct){
			correct=false; loop_num++;
			for(int i=0;i<edgesAdj.size();i++){
				curEdgeCate=edgesCate[i];
				if(curEdgeCate!=4 && edgesAdj.get(i).size() >= 2){// get a classified edge with at least two neighbors
					
					num_same=0; multiAdjCates=false; 
					for(int j=0; j<edgesAdj.get(i).size() && !multiAdjCates; j++){//check adjacent edges for current edge
						curEdgeAdj=edgesAdj.get(i).get(j);
						adjCate=edgesCate[curEdgeAdj];
						
						if (adjCate == curEdgeCate){
							num_same++;
						}
						else if(adjCate != 4){//find a classified neighboring edge with different category
							num_support=1; 
							j++;
							
							while(!multiAdjCates && j<edgesAdj.get(i).size()){
								curEdgeAdj=edgesAdj.get(i).get(j);
								if(edgesCate[curEdgeAdj] != 4){
									if(adjCate != edgesCate[curEdgeAdj]){
										multiAdjCates=true;
									}
									else{
										num_support++;
									}
								}
								j++;
							}
							
							//correct current edge category to be same as its classified neighbors
							if(!multiAdjCates && num_same==0 && num_support>=2){
								edgesCate[i] = adjCate;
								correct=true;
								correctEdges++;
							}
						}
					}
					
				}
			}
		}
		System.out.println("The correction step was repeated for "+loop_num+" times, and "+correctEdges+"  classified edges were corrected !");
//		System.out.println("finish correction steps !  ");
	}
	
    public static void propagation(ArrayList<ArrayList<Integer>> edgesAdj, int[] edgesCate){
		
		int curEdgeCate=-1, curEdgeAdj=-1, adjCate=-1, propagateEdges=0; 
		int num_support=0, num_same=0, loop_num=0;  
		boolean multiAdjCates=false, propagate=true; 
		
		while(propagate){
			propagate=false; loop_num++;
			for(int i=0;i<edgesAdj.size();i++){
				if(edgesCate[i]==4 && edgesAdj.get(i).size() >= 1){ //get an uncertain edge with at least one neighbor
					
					multiAdjCates=false; adjCate=-1; 
					for(int j=0; j<edgesAdj.get(i).size() && !multiAdjCates; j++){// check neighbors for current edge
						curEdgeAdj=edgesAdj.get(i).get(j);
						if(edgesCate[curEdgeAdj] != 4){//find a classified neighbors
							if(adjCate==-1){
								adjCate=edgesCate[curEdgeAdj]; 
							}
							else{
								if(adjCate != edgesCate[curEdgeAdj]){
									multiAdjCates=true;
								}
							}
						}												
					}
					
					if(!multiAdjCates && adjCate!=-1){ // propagate category of classified neighbors to the uncertain edge
						edgesCate[i] = adjCate;
						propagate=true;
						propagateEdges++;
					}
				}
			}
		}
		
		System.out.println("The propagation step was repeated for "+loop_num+" times, and "+propagateEdges+"  uncertain edges were classified !");
//		System.out.println("finish propagation steps !  ");
	}
    
    public static void finalContigCategoryFromEdges(String file_outdir, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges, int[] edgesLength, int[] edgesCate, boolean flyeAssembly){
    	try{
    		BufferedWriter writer = new BufferedWriter(new FileWriter(file_outdir+"3CAC_classification.fasta"));
    		int edge_ID=-1, edge_cate=-1, final_cate=-1; double sumLength=0, maxLength=0;
    		int[] cate_length=new int[4];
            String prefix="NODE_"; String finalContigCate="";
    		if(flyeAssembly){
    			prefix="contig_";
    		}
    		
    		for(int i=0; i<contigs.size(); i++){
    			for(int k=0;k<cate_length.length;k++){
    				cate_length[k]=0;
    			}
    			
    			maxLength=0; sumLength=0;
    			for(int j=0; j<contigEdges.get(i).size(); j++){
    				edge_ID=contigEdges.get(i).get(j)-1;
    				edge_cate=edgesCate[edge_ID]-1;
    				cate_length[edge_cate]=cate_length[edge_cate]+edgesLength[edge_ID];
    				sumLength=sumLength+edgesLength[edge_ID];
    				
    				if(maxLength<cate_length[edge_cate]){
    					maxLength=cate_length[edge_cate]; final_cate=edge_cate+1;
    				}
    			}
    			
    			//for contigs contain edges from different categories.
    			//if the categories with longest length doesn't account for 80% of the contig length, this contig will be classified as uncertain
    			if(maxLength/sumLength<0.8){ 
    				final_cate=4;
    			}
    			
    			if(final_cate==1){
    				finalContigCate="phage";
    			}
    			else if(final_cate==2){
    				finalContigCate="plasmid";
    			}
    			else if(final_cate==3){
    				finalContigCate="chromosome";
    			}
    			else{
    				finalContigCate="uncertain";
    			}
    			writer.write(prefix+contigs.get(i)+","+finalContigCate); writer.newLine();
    		}
    		writer.close();
    		System.out.println("The final classification of 3CAC can be found in: "+file_outdir+"3CAC_classification.fasta");
    		
    	}catch(Exception e){
			
		}
    }
		
	public static void main (String args[]){
		try{
			String assembler=""; // assembled by Flye or SPAdes
			String file_assembly_info=""; // file containing edges information for each contig
			String file_assembly_graph=""; // assembly graph file
			String file_initial=""; // initial classification file
			String file_outdir=""; // initial classification file
			
			for(int i=0; i<10; i=i+2){
				if(args[i].equalsIgnoreCase("--assembler")){
					assembler=args[i+1];
					System.out.println("assembler: "+assembler);
				}
				else if(args[i].equalsIgnoreCase("--initial")){
					file_initial=args[i+1];
					System.out.println("initial 3-way classification file: "+file_initial);
				}
				else if(args[i].equalsIgnoreCase("--graph")){
					file_assembly_graph=args[i+1];
					System.out.println("assembly graph file: "+file_assembly_graph);
				}
				else if(args[i].equalsIgnoreCase("--path")){
					file_assembly_info=args[i+1];
					System.out.println("contig path information file: "+file_assembly_info);
				}
				else if(args[i].equalsIgnoreCase("--output")){
					file_outdir=args[i+1];
					if(file_outdir.charAt(file_outdir.length()-1) != '/'){
						file_outdir=file_outdir+"/";
					}
					System.out.println("output directory: "+file_outdir);
				}
			}
			System.out.println();
			
			boolean flyeAssembly=false; 
			if(assembler.equalsIgnoreCase("flye") || assembler.equalsIgnoreCase("metaflye")){
				flyeAssembly=true;
			}
			
			
			//read edges in contigs
			int max_edge_ID=0;
			ArrayList<Integer> contigs=new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> contigEdges=new ArrayList<ArrayList<Integer>>();
			if(flyeAssembly){
				max_edge_ID = flyeContigsToEdgesInGraph(file_assembly_info, contigs, contigEdges);
			}
			else{
				max_edge_ID = spadesContigsToEdgesInGraph(file_assembly_info, contigs, contigEdges);
			}
			
			
			//read assembly graph to get adjacency and edge length
			int[] edgesLength=new int[max_edge_ID];
			int[] edgesCate=new int[max_edge_ID];
			ArrayList<ArrayList<Integer>> edgesAdj=new ArrayList<ArrayList<Integer>>();
			for(int i=0; i<edgesLength.length; i++){
				edgesLength[i]=-1; edgesCate[i]=4;
				ArrayList<Integer> tempAdj=new ArrayList<Integer>();
				edgesAdj.add(tempAdj);
			}
			readAssemblyGraph(file_assembly_graph, flyeAssembly, edgesLength, edgesAdj);
			
			
			//assign contig categories to edges 
			contigCategoryToEdges(file_initial, contigs, contigEdges, flyeAssembly, edgesCate);
			
			
			//correction and propagation step
			correction(edgesAdj, edgesCate);
			propagation(edgesAdj, edgesCate);
			
			
			//final contig category decided by edge category
			finalContigCategoryFromEdges(file_outdir, contigs, contigEdges, edgesLength, edgesCate, flyeAssembly);
			
			System.out.println("Finish 3CAC algorithm !");
    		
		}catch(Exception e){
			
		}
	}

}

