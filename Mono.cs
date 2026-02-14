/*
**  File: Mono.cs
**  Started: 
**  Contributors: Joanna Slusky, Meghan Franklin, Ryan Feehan et al.
**  Overview: 
**
**  About: 
**
**  Last Edited: 
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using betaBarrelProgram.BarrelStructures;
using System.Collections;
using System.IO;

namespace betaBarrelProgram
{
    namespace Mono
    {
        public class MonoProtein : Protein
        {
            public List<Chain> Chains { get; set; }
            /* public List<Chain> Chains = new List<BarrelStructures.Chain>(); */
            public int ChainCount { get; set; }
            public int totalResNum { get; set; }
            public string PdbName { get; set; }


            public MonoProtein(ref AtomParser.AtomCategory _myAtomCat, string PdbName)
            {
                /*this.Chains = new List<BarrelStructures.Chain>();*/
                this.ChainCount = 0;
                this.totalResNum = 0;
                this.PdbName = PdbName;
                this.Chains = new List<BarrelStructures.Chain>();

                for (int chainNum = 0; chainNum < _myAtomCat.chainAtomsList.Count; chainNum++)
                {
                    bool IsItProtein = false;
                    for (int atomNum = 0; atomNum < _myAtomCat.ChainAtomList[chainNum].cartnAtoms.Count; atomNum++)
                    {
                        if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomNum].atomName == "CA") IsItProtein = true;
                    }
                    if (IsItProtein == true)
                    {
                        Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, true, Global.DB_DIR, false);
                        //Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, false, Global.MONO_DB_DIR); //changed mono status for SS for solubles
                        this.Chains.Add(myChain);
                        this.ChainCount++;
                    }
                }
            }

            public MonoProtein(ref AtomParser.AtomCategory _myAtomCat, int chainNum, string PdbName)
            {
                this.ChainCount = 0;
                this.totalResNum = 0;
                this.PdbName = PdbName;
                this.Chains = new List<BarrelStructures.Chain>();
                BarrelStructures.Chain myChain = new BarrelStructures.Chain(ref _myAtomCat, chainNum, PdbName, true, Global.DB_DIR, false);
                this.Chains.Add(myChain);
            }

            public void translate(Vector3 translationCoords)
            {
                for (int chainCtr = 0; chainCtr < this.Chains.Count; chainCtr++)
                {
                    this.Chains[chainCtr].translate(translationCoords);
                }
            }
        }

        public class MonoBarrel : Barrel

        {
            public List<Strand> Strands { get; set; }
            public List<Vector3> NellipseCoords { get; set; }
            public List<Vector3> CellipseCoords { get; set; }
            public Vector3 Ncentroid { get; set; }
            public Vector3 Ccentroid { get; set; }
            public Vector3 Axis { get; set; }
            public bool Direction { get; set; }
            public double AvgTilt { get; set; }
            public double AvgTilt_even { get; set; }
            public double AvgTilt_odd { get; set; }
            public List<List<int>> protoBarrel { get; set; }
            public double AvgRadius { get; set; }
            public List<Res> LoopResies { get; set; }
            public string PdbName { get; set; }
            public string ChainName { get; set; }
            public List<double> StrandLength { get; set; }
            public Vector3 OriginalNcentroid { get; set; }
            public Vector3 OriginalCcentroid { get; set; }
            public Vector3 AxisVector { get; set; }
            public Vector3 NewCaxisPt { get; set; }
            public Vector3 OldCaxisPt { get; set; }
            public int ShearNum { get; set; }
            public List<double> PrevTwists { get; set; }
            public bool Success { get; set; }
            //public static string path = Global.OUTPUT_DIR;

            //barrel constructor 
            public MonoBarrel(Chain _myChain, Protein _myProtein)
            {
                this.Strands = new List<Strand>();
                this.ChainName = _myChain.ChainName;
                this.protoBarrel = new List<List<int>>();
                this.PdbName = _myChain.PdbName;
                this.Success = true;//need to actually check this at some point, but variable is currently only for all method
                string path = Global.OUTPUT_DIR;


                #region makeBarrel
                //Current method of defining strands
                createAllStrands(ref _myProtein);

                Console.WriteLine("After createStrands: \n");
                for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }
                //If DSSP definitions are important
                foreach (Res Res1 in _myChain){
                    addNeighs(Res1, ref _myChain);//just for adding neighbors right now (6/13/23)
				}
                concatinateStrands(ref _myChain);
                Console.WriteLine("After concat: \n");
                for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }
                makeBarrelCircular_2strandNeigh(ref _myChain);// keeps out the beta strands that occlude the barrel
                Console.WriteLine("After makeBarrelCircular(1st): \n");
                removeNonBarrelRes(ref _myChain);
                
             
                Console.WriteLine("After getting rid of loops: \n");
                for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }

                elongateStrands(ref _myChain);

                Console.WriteLine("After elongate strands: \n");
                for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }

                makeBarrelCircular_2strandNeigh(ref _myChain);

                //checkStrands is written to use protobarrel list before the strand list is created. -- has been deleted
                //redefine axis
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);
                this.Axis = this.Ncentroid - this.Ccentroid;
                /*end of 6H3I_A commented out*/

                for (int strandCtr = 0; strandCtr < protoBarrel.Count; strandCtr++)
                {
                    Strand newStrand = new Strand(_myChain, protoBarrel[strandCtr][0], protoBarrel[strandCtr][protoBarrel[strandCtr].Count - 1], strandCtr);
                    this.Strands.Add(newStrand);
                }

                #endregion
                SharedFunctions.writePymolScriptForStrands(this.Strands, Global.OUTPUT_DIR, Global.DB_DIR, this.PdbName);// added to look at SprA

                #region Amino Acid Pair
                try
                {
                    SharedFunctions.AminoAcidPairs(this.Strands, Global.OUTPUT_DIR, Global.DB_DIR, this.PdbName);
                }
                catch (Exception exception)
                {
                    Console.WriteLine("Failed to make pairs");
                    string fileLocation = Global.OUTPUT_DIR + "AAPairErrorLog.txt";
                    using (StreamWriter file = new StreamWriter(fileLocation, append: true))
                    {
                        string newLine = this.PdbName + "\t" + exception.InnerException + exception.Message + "\t" + exception.StackTrace + "\t" + exception.Source + "\t" + exception.Data + "\n";
                        file.WriteLine(newLine);
                    }
                }
                #endregion


                #region Rotate barrel to z-axis
                this.AvgRadius = SharedFunctions.setRadius(this.Strands, this.Axis, this.Ccentroid, this.Ncentroid);
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                this.AxisVector = SharedFunctions.getNormal(this.NellipseCoords, this.Ncentroid);

                this.OriginalCcentroid = this.Ccentroid;
                this.OriginalNcentroid = this.Ncentroid;
                this.NewCaxisPt = ((this.Ncentroid - this.Ccentroid).Length() * this.AxisVector) + this.Ncentroid;
                int axisDirection = 1;
                if ((this.NewCaxisPt - this.OriginalCcentroid).Length() > ((((this.Ncentroid - this.Ccentroid).Length() * -1 * this.AxisVector) + this.Ncentroid) - this.OriginalCcentroid).Length())
                {
                    axisDirection = -1;
                    this.NewCaxisPt = (((this.Ncentroid - this.Ccentroid).Length() * -1 * this.AxisVector) + this.Ncentroid);
                }

                this.Axis = this.NewCaxisPt - this.Ncentroid;
                this.OldCaxisPt = this.NewCaxisPt;
                this.rotateToZ(ref _myChain);

                if (this.Axis.Length() > 22.5) setBottom(ref _myChain);
                else centerZ(ref _myChain);

                rotate180(ref _myChain);

                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                this.AxisVector = SharedFunctions.getNormal(this.NellipseCoords, this.Ncentroid);
                this.NewCaxisPt = (axisDirection * (this.Ncentroid - this.Ccentroid).Length() * this.AxisVector) + this.Ncentroid;
                this.Axis = this.NewCaxisPt - this.Ncentroid;

                for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                {
                    for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count; resCtr++)
                    {
                        this.Strands[strandCtr].Residues[resCtr].Z = this.Strands[strandCtr].Residues[resCtr].BackboneCoords["CA"].Z;
                    }
                }
                #endregion

                SharedFunctions.setInOut(this.Strands, path, this.PdbName, this.Axis, this.Ccentroid, this.Ncentroid);

                /* ---------output information--------- */
	            this.StrandLength = SharedFunctions.getStrandLengths(this.Strands, path, this.PdbName);
                this.PrevTwists = SharedFunctions.writeTwists(this.Strands, Global.MONO_OUTPUT_DIR, this.PdbName);

                string use_dir = path + "ZCoords";
                if (!System.IO.Directory.Exists(use_dir))
                {
                    System.IO.Directory.CreateDirectory(use_dir);
                }

                string fileLocation15 = path + "ZCoords/AllZCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation15))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Strand" + "\t" + "Z-coord";
	                file.WriteLine(newLine);
	                foreach (Strand strand in this.Strands)
	                {
	                    foreach (Res res in strand)
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            if (atom.AtomName == "CA") file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", res.ThreeLetCode, res.SeqID, strand.StrandNum, strand.ChainName, atom.Coords.X, atom.Coords.Y, atom.Coords.Z);
	                        }
	                    }
	                }
	            }

	            //DSSP has to be used to define strands before the strands are created;

	            SharedFunctions.setInOut(this.Strands, path, this.PdbName, this.Axis, this.Ccentroid, this.Ncentroid);

	            this.AvgTilt = SharedFunctions.getTiltsByAA(this.Strands, path, this.PdbName, this.Axis, ref Global.AADict);
                this.Success = true;
                //getTyrVector(path, ref _myChain);
                WritePDB(ref _myChain);
                WriteDesignPos(ref _myChain);

            }

            public void createAllStrands(ref Protein _myProtein) //This does NOT create the list of Strands
            {
                bool strandStart = false;
                int strandStartRes = 0;
                int strandEndRes = 0;
                int CurrentResNum;
                int strandNum = 0;
                int betaScore = 0;
                List<List<int>> StrandsInChain = new List<List<int>>();
                List<int> myStrand = new List<int>();

                foreach (Chain chain in _myProtein.Chains)
                {
                    for (CurrentResNum = 0; CurrentResNum < chain.ResidueCount; CurrentResNum++)
                    {
                        if (strandStart == false)
                        {
                            if (CurrentResNum >= chain.ResidueCount - 5) continue;

                            if (chain.Residues[CurrentResNum].SSType == "B")
                            {
                                strandStartRes = chain.Residues[CurrentResNum].ResNum; //ResNum is index of residue w/in protein; Seq_ID is the PDB residue number
                                betaScore = 0;
                                for (int i = 1; i < 5; i++)
                                {
                                    if (chain.Residues[CurrentResNum + i].SSType == "B") betaScore++;
                                }
                                if (betaScore >= 2) strandStart = true;
                            }
                        }

                        else //Strand has been started
                        {
                            if ((CurrentResNum == chain.ResidueCount - 1) || (chain.Residues[CurrentResNum + 1].SeqID > (chain.Residues[CurrentResNum].SeqID + 7))) //If we are still forming a chain and reach the end of the residues OR if big gap in seqIDs, end the chain
                            {
                                goto EndChain;
                            }

                            else
                            {
                                if (chain.Residues[CurrentResNum].SSType == "B") continue;
                                else //if the chain is started and next residue is NOT beta-residue, determine if @ the end of a chain
                                {
                                    //if there's a big bend in chain, definitely end the chain
                                    if (SharedFunctions.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction) >= 80)
                                    {
                                       goto EndChain;
                                    }
                                    //if no bend in chain, check next residue in line 
                                    else
                                    {
                                        if (chain.Residues[CurrentResNum + 1].SSType == "B" && chain.Residues[CurrentResNum + 2].SSType == "B") continue; // XBB
                                        else if (chain.Residues[CurrentResNum + 1].SSType == "B" && chain.Residues[CurrentResNum + 2].SSType != "B") //XBX
                                        {
                                            double Angle1 = SharedFunctions.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
                                            double Angle2 = SharedFunctions.AngleBetween(chain.Residues[CurrentResNum + 1].Direction, chain.Residues[CurrentResNum + 2].Direction);
                                            //if (Angle1 > 71 || Angle2 > 71 || Angle1 < 45 || Angle2 < 45) goto EndChain; //B-type chain of 7ahl has angle at 71.5
                                            if (Angle1 > 71 || Angle2 > 71) goto EndChain; //B-type chain of 7ahl has angle at 71.5
                                            else continue;
                                        }
                                        else //XX
                                        {
                                            betaScore = 0;
                                            for (int i = 2; i < 5; i++)
                                            {
                                                try
                                                {
                                                    if (chain.Residues[CurrentResNum + i].SSType == "B") betaScore++;
                                                }
                                                catch (ArgumentOutOfRangeException)
                                                {
                                                    betaScore = 0;
                                                }

                                            }
                                            if (betaScore >= 3)
                                            {
                                                double Angle1 = SharedFunctions.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction);
                                                double Angle2 = SharedFunctions.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
                                                if (Angle1 + Angle2 > 129 || Angle1 > 80 || Angle2 > 80) //if there's a decent bend at the two non-beta residues between long beta-conf sequences //changed from 135 for 7AHL and 3W9T on 12/3 MWF
                                                {
                                                    goto EndChain;
                                                }/*
                                                double Angle3 = SharedFunctions.AngleBetween(chain.Residues[CurrentResNum-1].Direction, chain.Residues[CurrentResNum + 5].Direction);
                                                if (Angle3 > 140) //if there's a decent bend at the two non-beta residues between long beta-conf sequences //changed from 135 for 7AHL and 3W9T on 12/3 MWF
                                                {
                                                    goto EndChain;
                                                }
                                                */
                                            }
                                            else //this means two res in a row were non-beta AND that less than three res of next 3 were also not beta 
                                            {
                                                goto EndChain;
                                            }

                                        }
                                    }
                                } //end of dealing with non-Beta conformation residues
                            }
                        } //end of strandStart == true
                        continue;
                    EndChain:
                        {
                            strandEndRes = CurrentResNum - 1;
                            if ((strandEndRes - strandStartRes) >= 2)
                            {
                                for (int j = strandStartRes; j <= strandEndRes; j++)
                                {
                                    myStrand.Add(chain.Residues[j].ResNum);
                                    //chain.Residues[j].SSType = "B";
                                }
                                List<int> newList = new List<int>();
                                newList.AddRange(myStrand);
                                StrandsInChain.Add(newList);
                                myStrand.Clear();
                                strandNum++;
                            }
                            strandStart = false;
                        }

                    }

                    foreach (List<int> this_strand in StrandsInChain){protoBarrel.Add(this_strand); }
                    StrandsInChain.Clear();
                    strandNum = 0;

                }

            //for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myProtein.Chains[0].Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myProtein.Chains[0].Residues[this.protoBarrel[x].Last()].SeqID); }

            }//end of create strands fxn

            public void concatinateStrands(ref Chain _myChain)
            {
                int nextStrand = 0;
                for (int strandNum = 0; strandNum < this.protoBarrel.Count - 1; strandNum++)
                {
                    nextStrand = strandNum + 1;
                    for (int resNum1 = 0; resNum1 < this.protoBarrel[strandNum].Count; resNum1++)
                    {
                        for (int nCtr = 0; nCtr < _myChain.Residues[this.protoBarrel[strandNum][resNum1]].Neighbors.Count; nCtr++) //goes through neighbors of every residue in the current strand
                        {
                            for (int findStrandNeigh = 0; findStrandNeigh < this.protoBarrel.Count - 1; findStrandNeigh++)
                            {
                                if (findStrandNeigh != strandNum && findStrandNeigh != nextStrand && this.protoBarrel[findStrandNeigh].Contains(_myChain.Residues[this.protoBarrel[strandNum][resNum1]].Neighbors[nCtr])) //finds a neighboring strand to current strand
                                {
                                    for (int nextResNum = 0; nextResNum < this.protoBarrel[nextStrand].Count; nextResNum++)
                                    {
                                        for (int nextNCtr = 0; nextNCtr < _myChain.Residues[this.protoBarrel[nextStrand][nextResNum]].Neighbors.Count; nextNCtr++)
                                        {
                                            if (this.protoBarrel[findStrandNeigh].Contains(_myChain.Residues[this.protoBarrel[nextStrand][nextResNum]].Neighbors[nextNCtr])) // if the neighboring strand of current strand contains neighboring residues with next strand, concatinate current and next strand 
                                            {
                                                Vector3 strand1vec = _myChain.Residues[protoBarrel[strandNum].Last()].BackboneCoords["CA"] - _myChain.Residues[protoBarrel[strandNum][0]].BackboneCoords["CA"];
                                                Vector3 strand2vec = _myChain.Residues[protoBarrel[nextStrand].Last()].BackboneCoords["CA"] - _myChain.Residues[protoBarrel[nextStrand][0]].BackboneCoords["CA"];
                                                double angle = SharedFunctions.AngleBetween(strand1vec, strand2vec);
                                                if (angle <= 80)
                                                {
                                                    protoBarrel[strandNum].AddRange(protoBarrel[nextStrand]);
                                                    protoBarrel.RemoveRange(nextStrand, 1);
                                                    Console.WriteLine("Combined two strands");

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            public void elongateStrands(ref Chain _myChain)
            {
                for (int strandNum = 0; strandNum < this.protoBarrel.Count - 1; strandNum++)
                {
                    int firstRes = this.protoBarrel[strandNum][0];
                    while (_myChain.Residues[firstRes - 1].SSType == "B") //adds previous residues with beta conformation
                    {
                        if (this.protoBarrel[strandNum].Contains(firstRes - 1) == false)
                        {
                            protoBarrel[strandNum].Add(_myChain.Residues[firstRes - 1].ResNum);
                            protoBarrel[strandNum].Sort();
                        }
                        firstRes--;
                    }
                    int lastRes = this.protoBarrel[strandNum][this.protoBarrel[strandNum].Count - 1];
                    while (_myChain.Residues[lastRes + 1].SSType == "B") //adds ending residues with beta confirmation
                    {
                        if (this.protoBarrel[strandNum].Contains(lastRes + 1) == false)
                        {
                            protoBarrel[strandNum].Add(_myChain.Residues[lastRes + 1].ResNum);
                            protoBarrel[strandNum].Sort();
                        }
                        lastRes++;
                    }
                }

            }

            public void checkStrandDefnsDSSP(ref Chain _myChain) //Added 6-5-17 for loops
	        {
	            for (int strandNum = 0; strandNum < this.protoBarrel.Count; strandNum ++)
	            {
	                //Remove turn residues
	                while (_myChain.Residues[this.protoBarrel[strandNum][0]].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum][0]].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(0);
	                    Console.WriteLine("Removed res{1} from beg of strand {0}", strandNum, removed_res);
	                }

	                while (_myChain.Residues[this.protoBarrel[strandNum].Last()].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum].Last()].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(this.protoBarrel[strandNum].Count - 1);
	                    Console.WriteLine("Removed res{1} from end of strand {0}", strandNum, removed_res);
	                }


	                //Add to the beginning of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum][0] - 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum][0] - 1 > this.protoBarrel[strandNum - 1].Last() + 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum][0] - 1], ref _myChain);
	                                Console.WriteLine("added res{1} to beg of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break;  }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into prev strand"); break; }
                        
	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }

	                //For adding to the ends of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum].Last() + 1 < this.protoBarrel[strandNum + 1][0] - 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1], ref _myChain);
	                                Console.WriteLine("added res{1} to end of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break;  }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into next strand"); break; }
	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }
	            }

	        }

            public void addNeighs(Res Res1, ref Chain _myChain) //Added 6-5-17 for loops
            {
                double minD = 3.0;

                foreach (Res Res2 in _myChain)
                {
                    if (Math.Abs(Res1.ResNum - Res2.ResNum) > 3)
                    {
                        double d = (Res1.Atoms[0].Coords - Res2.BackboneCoords["O"]).Length();
                        if (d < minD)// && (Res2.SSType == "B" || Res2.DSSP == "E"))
                        {
                            if (Res1.Neighbors.Contains(Res2.ResNum) == false) Res1.Neighbors.Add(Res2.ResNum);
                            if (Res2.Neighbors.Contains(Res1.ResNum) == false) Res2.Neighbors.Add(Res1.ResNum);
                        }
                    }
                }
            }

            public void checkDSSPNeighs(Res Res1, ref Chain _myChain) //Added 6-5-17 for loops
	        {
	            double minD = 2.75;

	            foreach (Res Res2 in _myChain)
	            {
	                if (Math.Abs(Res1.SeqID - Res2.SeqID) > 2)
	                {
	                    double d = (Res1.Atoms[0].Hydrogen - Res2.BackboneCoords["O"]).Length();
	                    if (d < minD && (Res2.SSType == "B" || Res2.DSSP == "E"))
	                    {
	                        if (Res1.Neighbors.Contains(Res2.ResNum) == false) Res1.Neighbors.Add(Res2.ResNum);
	                        if (Res2.Neighbors.Contains(Res1.ResNum) == false) Res2.Neighbors.Add(Res1.ResNum);
	                    }
	                }
	            }                       
	        }

            // this will calculate the shear number of the beta barrel.  see murzin lesk and chothia 1994 http://www.mrc-lmb.cam.ac.uk/tcb/pdf/chc/97_jmb_236_1369_94.pdf
            public int shearNum(ref Chain _myChain)
            {
                string filelocation = Global.OUTPUT_DIR + "ShearNum_" + PdbName + ".txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(filelocation))
                {
                    // first determine nearest hbonding neighbor with N coming from N-term  (arbitrary)           
                    for (int strandCtr = 0; strandCtr < this.Strands.Count(); strandCtr++)
                    {
                        int strand2 = strandCtr + 1;
                        if (strandCtr == this.Strands.Count() - 1) strand2 = 0;

                        for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count(); resCtr++)
                        {

                            Res firstRes = this.Strands[strandCtr].Residues[resCtr];
                            double minD = 3.4;
                            if (firstRes.Atoms[0].AtomName == "N" && firstRes.ThreeLetCode != "PRO")
                            {
                                Atom firstAtom = firstRes.Atoms[0];

                                for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
                                {
                                    Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[3];
                                    if (secondAtom.AtomName == "O")
                                    {
                                        double d = (secondAtom.Coords - firstAtom.Coords).Length();
                                        if (d < minD)
                                        {
                                            minD = d;
                                            firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
                                            firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
                                            firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length();
                                        }

                                    }

                                }


                            }
                            // looking at o=n as well just to be safe... 
                            if (firstRes.h_bonder == 0 && firstRes.Atoms[3].AtomName == "O")
                            {
                                Atom firstAtom = firstRes.Atoms[3];

                                for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
                                {
                                    Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[0];
                                    if (secondAtom.AtomName == "N")
                                    {
                                        double d = (secondAtom.Coords - firstAtom.Coords).Length();

                                        if (d < minD)
                                        {
                                            minD = d;
                                            firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
                                            firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
                                            firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length();
                                        }
                                    }
                                }
                            }
                            if (firstRes.h_bonder == 0) file.WriteLine("strand {0} res {1} has no partner!", strandCtr, firstRes.SeqID);
                            else file.WriteLine("strand {0} res {1}s partner is {2}; k = {3}", strandCtr, firstRes.SeqID, firstRes.h_bonderID, firstRes.ResStrandNum);
                        }
                    }
                    return 1;
                }
            }


            public void rotate180(ref Chain _myChain)
            {
                double phi = Math.PI;
                Matrix4x4 rotationMatrixY2 = new Matrix4x4((float)Math.Cos(phi), 0, (float)Math.Sin(phi), 0,
                                                         0, 1, 0, 0,
                                                           (float)(-1 * Math.Sin(phi)), 0, (float)Math.Cos(phi),
                                                         0, 0, 0, 0, 0);

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = Vector3.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords, rotationMatrixY2);
                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
            }

            public void setBottom(ref Chain _myChain)
            {
                double adjuster = 12.5;
                //double adjuster = .0962 * this.Axis.Length() + 8.3358;
                //double adjuster = 11.75;
                // if (this.Axis.Length() < 30) adjuster = 11;
                // if (this.Axis.Length() > 40) adjuster = 12.75;
                double adjustment = adjuster - this.Ncentroid.Z;

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {

                        Vector3 newCoords = new Vector3(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.X, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Y, (float)(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.Z + (adjustment)));
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = newCoords;


                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
            }

            public void centerZ(ref Chain _myChain)
            {
                double avgZ = (this.Ncentroid.Z + this.Ccentroid.Z) / 2;

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {

                        Vector3 newCoords = new Vector3(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.X, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Y, (float)(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.Z - (avgZ)));
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = newCoords;


                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }


            }

            public void rotateToZ(ref Chain _myChain)
            {
                double length = this.Axis.Length();

                double psi = Math.Atan(this.Axis.Y / this.Axis.X);
                if (0 == this.Axis.X) { psi = 0; }
                Matrix4x4 rotationMatrixZT = new Matrix4x4((float)Math.Cos(psi), (float)Math.Sin(psi), 0, 0, -1 * (float)Math.Sin(psi), (float)Math.Cos(psi), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
                Matrix4x4 rotationMatrixZ = new Matrix4x4((float)Math.Cos(psi), -1 * (float)Math.Sin(psi), 0, 0, (float)Math.Sin(psi), (float)Math.Cos(psi), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
                Vector3 axis1 = new Vector3();
                axis1 = Vector3.Transform(this.Axis, rotationMatrixZ);

                double theta = Math.Atan(axis1.X / axis1.Z);
                //if (0 == this.Axis.Z) { theta = 0; }
                Matrix4x4 rotationMatrixYT = new Matrix4x4((float)Math.Cos(theta), 0, -1 * (float)Math.Sin(theta), 0, 0, 1, 0, 0, (float)Math.Sin(theta), 0, (float)Math.Cos(theta), 0, 0, 0, 0, 0);
                Matrix4x4 rotationMatrixY = new Matrix4x4((float)Math.Cos(theta), 0, (float)Math.Sin(theta), 0, 0, 1, 0, 0, -1 * (float)Math.Sin(theta), 0, (float)Math.Cos(theta), 0, 0, 0, 0, 0);

                Vector3 axis2 = new Vector3();
                axis2 = Vector3.Transform(axis1, rotationMatrixY);


                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = Vector3.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords, rotationMatrixZ);
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = Vector3.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords, rotationMatrixY);



                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }

                if (this.Strands[1].Residues[0].BackboneCoords["CA"].Z > this.Strands[0].Residues[0].BackboneCoords["CA"].Z)
                {
                    double phi = Math.PI;
                    Matrix4x4 rotationMatrixY2 = new Matrix4x4((float)Math.Cos(phi), 0, (float)Math.Sin(phi), 0, 0, 1, 0, 0, -1 * (float)Math.Sin(phi), 0, (float)Math.Cos(phi), 0, 0, 0, 0, 0);
                    Vector3 axis3 = new Vector3();
                    axis3 = Vector3.Transform(axis2, rotationMatrixY2);


                    for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                    {
                        for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                        {
                            _myChain.Residues[resCtr].Atoms[atomCtr].Coords = Vector3.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords, rotationMatrixY2);
                            if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                            {

                                _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                            }
                        }
                    }
                }
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                Vector3 translate = new Vector3(this.Ccentroid.X, this.Ccentroid.Y, 0);

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords -= translate;

                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);
            }

            public void makeBarrelCircular_2strandNeigh(ref Chain _myChain)
            {

                for (int strndCtr2 = 0; strndCtr2 < this.protoBarrel.Count; strndCtr2++)
                {
                    if (this.protoBarrel[strndCtr2].Count == 0) this.protoBarrel.RemoveAt(strndCtr2);
                    for (int i = this.protoBarrel[strndCtr2][0]; i < this.protoBarrel[strndCtr2][this.protoBarrel[strndCtr2].Count - 1]; i++)
                    {
                        if (this.protoBarrel[strndCtr2].Contains(i) == false)
                        {
                            this.protoBarrel[strndCtr2].Add(i);
                            this.protoBarrel[strndCtr2].Sort();
                        }


                    }
                }

                bool allStrands2Neighs = true;
                do
                {
                    allStrands2Neighs = true;
                    for (int strndCtr1 = 0; strndCtr1 < this.protoBarrel.Count;)
                    {
                        List<int> myStrandNeighList = new List<int>();

                        for (int strndCtr2 = 0; strndCtr2 < this.protoBarrel.Count; strndCtr2++)
                        {
                            if (strndCtr1 != strndCtr2)
                            {
                                for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr1].Count; resCtr++)
                                { // for each residue
                                    for (int nCtr = 0; nCtr < _myChain.Residues[this.protoBarrel[strndCtr1][resCtr]].Neighbors.Count; nCtr++)
                                    {
                                        if (this.protoBarrel[strndCtr2].Contains(_myChain.Residues[this.protoBarrel[strndCtr1][resCtr]].Neighbors[nCtr]))
                                        {
                                            if (myStrandNeighList.Contains(strndCtr2) == false) myStrandNeighList.Add(strndCtr2);
                                        }
                                    }
                                }
                            }

                        }
                        if (myStrandNeighList.Count< 2)
                        {
                            Console.WriteLine("removing strand {0}", strndCtr1);
                            this.protoBarrel.Remove(this.protoBarrel[strndCtr1]);
                            allStrands2Neighs = false;
                        }
                        else strndCtr1++;
                    }


                    } while (false== allStrands2Neighs);

                }

            public void removeNonBarrelRes(ref Chain _myChain)
            {// removes residues that are not h-bonded to the next or the previous strand.
             //bool deletedSomething = false;
                for (int strndCtr = 0; strndCtr < this.protoBarrel.Count; strndCtr++)
                {
                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count;)
                    {
                        bool markForDeletion = true;

                        if (_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors.Count != 0)

                        {
                            for (int nCtr = 0; nCtr < _myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors.Count; nCtr++)
                            {
                                if (strndCtr == 0)
                                {
                                    if (this.protoBarrel[this.protoBarrel.Count - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else if (strndCtr == this.protoBarrel.Count - 1)
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[0].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                            }
                        }

                        if (markForDeletion == true)
                        {
                            this.protoBarrel[strndCtr].RemoveAt(resCtr);
                        }
                        else resCtr++;
                        //resCtr--;

                    }

                    this.protoBarrel[strndCtr].Sort();

                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count;)
                    {// if there are no residues before or after 5 that are part of the sheet delete that residue too.
                        bool markForDeletion = true;
                        for (int ctr2 = this.protoBarrel[strndCtr][resCtr] - 5; ctr2 < this.protoBarrel[strndCtr][resCtr] + 5; ctr2++)
                        {
                            if (ctr2 > -1 && ctr2 != this.protoBarrel[strndCtr][resCtr])
                            {
                                if (this.protoBarrel[strndCtr].Contains(ctr2) == true)
                                {
                                    markForDeletion = false;
                                }

                            }
                        }
                        if (markForDeletion == true)
                        {
                            this.protoBarrel[strndCtr].RemoveAt(resCtr);
                            Console.WriteLine("removed a residue from Strand {0} because no nearby residues", strndCtr);
                            //deletedSomething = true;
                        }
                        else resCtr++;
                    }

                    if (this.protoBarrel[strndCtr].Count == 0) protoBarrel.RemoveAt(strndCtr);
                }
            }

            public void setCEllipseCoords(ref Chain _myChain)
            {
                //this is the top (extracellular) ellipse
                List<Vector3> myEllipse = new List<Vector3>();
                Vector3 centroid = new Vector3();

                for (int strandCtr = 0; strandCtr < this.protoBarrel.Count; strandCtr++)
                {
                    Vector3 firstCA = new Vector3();
                    firstCA = _myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"];
                    Vector3 lastCA = new Vector3();
                    lastCA = _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"];



                    if (strandCtr % 2 == 0)
                    {
                        myEllipse.Add(lastCA);
                        centroid += lastCA;

                    }
                    else
                    {
                        myEllipse.Add(firstCA);
                        centroid += firstCA;

                    }

                }
                this.CellipseCoords = myEllipse;
                centroid = centroid / this.protoBarrel.Count;
                this.Ccentroid = centroid;


            }

            public void setNEllipseCoords(ref Chain _myChain)
            {
                //this is the bottom (periplasmic) ellipse
                List<Vector3> myEllipse = new List<Vector3>();
                Vector3 centroid = new Vector3();

                for (int strandCtr = 0; strandCtr < this.protoBarrel.Count; strandCtr++)
                {
                    Vector3 firstCA = new Vector3();
                    firstCA = _myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"];
                    Vector3 lastCA = new Vector3();
                    lastCA = _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"];


                    if (strandCtr % 2 == 0)
                    {
                        myEllipse.Add(firstCA);
                        centroid += firstCA;

                    }
                    else
                    {
                        myEllipse.Add(lastCA);
                        centroid += lastCA;
                    }

                }
                this.NellipseCoords = myEllipse;
                centroid = centroid / this.protoBarrel.Count;
                this.Ncentroid = centroid;


            }

            public void WritePDB(ref Chain _myChain)
            {
                string fileName = this.PdbName +"_"+this.ChainName+ ".pdb";
                string outputDir = Global.OUTPUT_DIR + "axisAlignedPDBS/";
                SharedFunctions.create_dir(outputDir);
                FileStream fileStream = new FileStream(Path.Combine(outputDir, fileName), FileMode.Create, FileAccess.Write);

                StreamWriter fileWriter = new StreamWriter(fileStream);
                string header = "HEADER    " + this.PdbName + "_" + this.ChainName + "                    " + DateTime.Now;
                fileWriter.WriteLine(header);

                try
                {
                    string line = "";
                    int atomCount = 1;

                    
                        foreach (Res res in _myChain)
                        {
                            foreach (Atom atom in res.Atoms)
                            {
                                line = "ATOM  ";

                                string atomIdStr = atomCount.ToString();
                                line += atomIdStr.PadLeft(5, ' ');
                                line += " ";

                                string atomName = atom.AtomName;
                                if (atomName != "" && atom.AtomType != "H" && atomName.Length < 4)
                                {
                                    atomName = " " + atomName;
                                }
                                line += atomName.PadRight(4, ' ');

                                line += " ";
                                line += res.ThreeLetCode;

                                line += " ";
                                line += res.ChainName;

                                line += res.SeqID.ToString().PadLeft(4, ' ');
                                line += "    ";
                                line += FormatDoubleString(atom.Coords.X, 4, 3);
                                line += FormatDoubleString(atom.Coords.Y, 4, 3);
                                line += FormatDoubleString(atom.Coords.Z, 4, 3);

                                line += "  1.00"; //(for dummy occupancy)

                                line += "  0.00"; //(for dummy bfactor)

                                line += "    ";
                                line += atom.AtomType;
                                fileWriter.WriteLine(line);
                                atomCount++;
                            }
                        }
                        fileWriter.WriteLine("END");
                 
                }
                catch (Exception ex)
                {
                    string errorMsg = ex.Message;
                    throw ex;
                }
                finally
                {
                    fileWriter.Close();
                }
            }

            private string FormatDoubleString(float val, int numPre, int numPost)
            {
                string valStr = val.ToString();
                int dotIndex = valStr.IndexOf(".");
                if (dotIndex == -1)
                {
                    // return the int part, plus ".0  "
                    valStr = valStr.PadLeft(numPre, ' ');
                    valStr += ".";
                    int i = 0;
                    while (i < numPost)
                    {
                        valStr += "0";
                        i++;
                    }
                    return valStr;
                }
                string intPartStr = valStr.Substring(0, dotIndex).PadLeft(numPre, ' ');
                int subStrLen = valStr.Length - dotIndex - 1;
                if (subStrLen > numPost)
                {
                    subStrLen = numPost;
                }
                string fractStr = valStr.Substring(dotIndex + 1, subStrLen).PadRight(3, '0');
                return intPartStr + "." + fractStr;
            }

            public IEnumerator<Strand> GetEnumerator()
            {
                return Strands.GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return this.GetEnumerator();
            }

            public void WriteDesignPos(ref Chain _myChain)
            {
                string fileName = this.PdbName.ToLower() + "_" + this.ChainName + ".pos";
                string outputDir = Global.OUTPUT_DIR + "DesignPos/";
                SharedFunctions.create_dir(outputDir);
                FileStream fileStream = new FileStream(Path.Combine(outputDir, fileName), FileMode.Create, FileAccess.Write);

                StreamWriter fileWriter = new StreamWriter(fileStream);

                foreach (Res res in _myChain)
                {
                    if (res.BackboneCoords["CA"].Z >= -1)
                    {
                        //resNum because Rosetta will renumber
                        fileWriter.Write(res.ResNum.ToString() + " ");
                    }
                }
                fileWriter.Close();
            }
        }
	
    }

 }