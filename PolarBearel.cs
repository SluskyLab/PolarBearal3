

using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Linq;
using betaBarrelProgram.BarrelStructures;

namespace betaBarrelProgram
{

    public class PolarBearal
    {
        //polarBearal Variable 
        public static string PolarBearal_OUTPUT_DIR = Global.OUTPUT_DIR + "PolarBearal/";
        public static string PolarBearal_INPUT_DB_FILE = Global.DB_file; //Global.POLARBEARAL_DIR + "/DB/MonoDB_v5.txt"; //Global.MONO_DB_file;

        double zone = 13;// +-zone (a.k.a. membrane) is considered when use_zone==true

        static public void menu()
        {
            Console.WriteLine("1. Generate empty result files");
            Console.WriteLine("2. Create Protein Database");
            Console.WriteLine("3. Run Protein Database");
            Console.WriteLine("4. ?");
            Console.WriteLine("10. Quit");
        }

        static public void RunPolarBearal()
        {
            PolarBearal_OUTPUT_DIR = Global.OUTPUT_DIR + "PolarBearal/";
            PolarBearal_INPUT_DB_FILE = Global.DB_file;

            string choice = "";
            while (choice != "10")
            {
                menu();
                choice = Console.ReadLine();
                switch (choice)
                {
                    case "1":
                        PrepResults();
                        break;
                    case "2":
                        using (StreamWriter log = File.AppendText(Global.OUTPUT_DIR + "log.txt"))
                        {
                            log.WriteLine("starting create betabarrel ProteinDB at {0}", DateTime.Now);
                            log.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", "PdbName", "Success", "barrelType", "chainID", "numberOfStrands");
                        }
                        CreateBetaBarrelProteinDatabase();
                        using (StreamWriter log = File.AppendText(Global.OUTPUT_DIR + "log.txt"))
                        {
                            log.WriteLine("ending create betabarrel ProteinDB at {0}", DateTime.Now);
                        }
                        break;
                    case "3":
                        RunBetaBarrelProteinDatabase();
                        Results();
                        break;
                    case "4":
                        break;
                    default:
                        choice = "10";
                        break;
                }
            }
        }

        static public void PrepResults()
        {
            string output_sub_dir;

            Console.WriteLine("Currently in {0}", System.IO.Directory.GetCurrentDirectory());
            Console.WriteLine("Output will be in {0}", PolarBearal_OUTPUT_DIR);

            Console.WriteLine("Would you like to delete all files currently in output dir (Y/N)?");
            string response = Console.ReadLine();
            if (response == "Y")
            {
                if (System.IO.Directory.Exists(PolarBearal_OUTPUT_DIR))
                {
                    System.IO.Directory.Delete(PolarBearal_OUTPUT_DIR, true);
                }
            }

            // build any needed paths for output
            if (!System.IO.Directory.Exists(PolarBearal_OUTPUT_DIR))
            {
                System.IO.Directory.CreateDirectory(PolarBearal_OUTPUT_DIR);
            }
            output_sub_dir = PolarBearal_OUTPUT_DIR + "betaBarrel_aaOnly";
            if (!System.IO.Directory.Exists(output_sub_dir))
            {
                System.IO.Directory.CreateDirectory(output_sub_dir);
            }
            output_sub_dir = PolarBearal_OUTPUT_DIR + "betaBarrelRawData";
            if (!System.IO.Directory.Exists(output_sub_dir))
            {
                System.IO.Directory.CreateDirectory(output_sub_dir);
            }
            output_sub_dir = PolarBearal_OUTPUT_DIR + "betaBarrelStrands";
            if (!System.IO.Directory.Exists(output_sub_dir))
            {
                System.IO.Directory.CreateDirectory(output_sub_dir);
            }
            output_sub_dir = PolarBearal_OUTPUT_DIR + "betaBarrelStrands_UseForCalc";
            if (!System.IO.Directory.Exists(output_sub_dir))
            {
                System.IO.Directory.CreateDirectory(output_sub_dir);
            }
            //create (or recreate) blank result files with headers
            using (StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + "PolarBearalResults.txt")) { }
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + "ExamineZ.txt")) { }
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + "PinNoutByProtein.txt"))
            {
                output.Write("\n {0} \t {1} \t {2} \t {3} \t {4} \t {5}", "PDB", "AAs", "IN", "Pin", "OUT", "Nout");
            }
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + "DisplayAngles.txt"))
            {
                output.Write("{0}\t{1}\t{2}", "aa", "angle", "inward facing");
            }
        }

       
        static public void CreateBetaBarrelProteinDatabase()
        {
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();

            string fileOfPDBs = PolarBearal_INPUT_DB_FILE;
            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    string fileLocation2 = PolarBearal_OUTPUT_DIR + "AllBarrelChar.txt";
                    using (System.IO.StreamWriter AllBarrel_output = new System.IO.StreamWriter(fileLocation2))
                    {
                        //string newLine = "PDB" + "\t\t" + "Total Strands" +"\t" + "Length" + "\t" + "AvgLength" + "\t" + "MinLength" + "\t" + "MaxLength" + "\t" + "Radius" + "\t" + "Barrel Tilt";
                        //file.WriteLine(newLine);
                        AllBarrel_output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", "PDB", "total_strands", "strand_cnt", "avg_length", "min_length", "max_length", "avg_radius", "avg_tilt", "shear_num", "start_res", "end_res");
                        // Read and display lines from the file until the end of the file is reached.
                        while ((line = sr.ReadLine()) != null)
                        {
                            string[] splitLine = line.Split(new char[] { ' ', '\t', ',' });
                            string pdb = splitLine[0];
                            //Console.Write(pdb);

                            //if (pdb != "IDs")
                            try
                            {
                                string fileName = pdb;
                                //string fileName = pdb + ".pdb";
                                //Barrel myBarrel = Program.runThisBetaBarrel(pdb, Global.METHOD);
                                BarrelStructures.Protein _protein = null;
                                BarrelStructures.Barrel myBarrel = null;
                                Program.runThisBetaBarrel(pdb, Global.METHOD, ref myBarrel, ref _protein);

                                SharedFunctions.LogBarrel(ref myBarrel, Global.METHOD);
                                PolarBearal roar = new PolarBearal(ref myBarrel);

                                string PDB = myBarrel.PdbName;
                                string total_strands = myBarrel.Axis.Length().ToString();
                                string avg_length = myBarrel.StrandLength.Average().ToString();
                                string min_length = myBarrel.StrandLength.Min().ToString();
                                string max_length = myBarrel.StrandLength.Max().ToString();
                                string avg_radius = myBarrel.AvgRadius.ToString();
                                string strand_cnt = myBarrel.Strands.Count.ToString();
                                string avg_tilt = myBarrel.AvgTilt.ToString();
                                string shear_num = myBarrel.ShearNum.ToString();
                                string start_res = myBarrel.Strands[0].ResNumStart.ToString();
                                string end_res = myBarrel.Strands[myBarrel.Strands.Count() - 1].ResNumEnd.ToString();
                                AllBarrel_output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", PDB, total_strands, strand_cnt, avg_length, min_length, max_length, avg_radius, avg_tilt, shear_num, start_res, end_res);


                            }
                            catch
                            {
                                AllBarrel_output.WriteLine("FAILED{0}", pdb);
                            }
                        }
                    }
                }
                //Console.WriteLine("Number of Proteins: {0} \t AAs: {1} \t Double Checked Directions: {2}", totalProteins, totalAAs, numDoubleChecks);
            }
            else
            {
                Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
                Console.WriteLine("could not open {0}", fileOfPDBs);
                Console.ReadLine();

            }
        }

        static public void RunBetaBarrelProteinDatabase()
        {
            string fileOfPDBs = PolarBearal_INPUT_DB_FILE;

            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        string[] splitLine = Array.FindAll<string>(((string)line).Split(
                            new char[] { ' ', '\t', ',' }), delegate (string s) { return !String.IsNullOrEmpty(s); });
                        string pdb = splitLine[0];
                        if (pdb != "IDs")
                        {
                            string fileName = pdb.ToUpper(); //RYANpdb.Substring(0, 4).ToUpper();

                            PolarBearal polarRetest = new PolarBearal(fileName);
                        }
                    }
                }
            }
            else
            {
                Console.WriteLine("could not open {0}", fileOfPDBs);
                Console.ReadLine();

            }
        }

        static public void Results()
        {
            string file = @"PolarBearalResults.txt";
            using (StreamWriter output = File.AppendText(PolarBearal_OUTPUT_DIR + file))
            {
                output.Write("\n\nProteins:\t{0} \t\t Strands:\t{1} \t\t AAs:\t{2} \t\t used AAs:\t{3}", PolarBearal._totalProteins, PolarBearal._totalStrands, PolarBearal._totalAAs, PolarBearal._usedAAs);
                output.Write("\n\n\n");

                int num = 0;
                output.Write("Strand Sizes:\n");
                output.Write("Length\t");
                foreach (int i in PolarBearal._strandSizes) output.Write("{0}\t", num++);
                output.Write("\nObserved\t");
                foreach (int i in PolarBearal._strandSizes) output.Write("{0}\t", i);
                output.Write("\n\n\n");

                
                num = 0;
                output.Write("Strand sizes within zone:\n");
                output.Write("Length\t");
                foreach (int i in PolarBearal._AlternatingStrandSeqLengths) output.Write("{0}\t", num++);
                output.Write("\nObserved\t");
                foreach (int i in PolarBearal._AlternatingStrandSeqLengths) output.Write("{0}\t", i);
                output.Write("\n\n\n");

                num = 0;
                output.Write("Alternating Direction Seq. Lengths:\n");
                output.Write("Length\t");
                foreach (int i in PolarBearal._AlternatingInOutSeqLengths) output.Write("{0}\t", num++);
                output.Write("\nObserved\t");
                foreach (int i in PolarBearal._AlternatingInOutSeqLengths) output.Write("{0}\t", i);
                output.Write("\n\n\n");

                num = 0;
                output.Write("Alternating Hydrophobicity Seq. Lengths:\n");
                output.Write("Length\t");
                foreach (int i in PolarBearal._AlternatingPolNonPSeqLengths) output.Write("{0}\t", num++);
                output.Write("\nObserved\t");
                foreach (int i in PolarBearal._AlternatingPolNonPSeqLengths) output.Write("{0}\t", i);
                output.Write("\n\n\n");
                num = 0;
                output.Write("Alternating Direction and Hydrophobicity Seq. Lengths:\n");
                output.Write("Length\t");
                foreach (int i in PolarBearal._AlternatingBothSeqLengths ) output.Write("{0}\t", num++);
                output.Write("\nObserved\t");
                foreach (int i in PolarBearal._AlternatingBothSeqLengths) output.Write("{0}\t", i);
                output.Write("\n\n\n");

                output.Write("AA directions\t");
                output.Write("\n\t\tL \tW \tP \tV \tI \tF \tY \tH \tA \tM \tG \tT \tQ \tD \tN \tK \tS \tR \tE \tC");
                output.Write("\nTotal\t");
                foreach (int i in __aa) output.Write("{0}\t", i);
                output.Write("\nInward\t");
                foreach (int i in _aaInward) output.Write("{0}\t", i);
                output.Write("\nOutward\t");
                foreach (int i in _aaOutward) output.Write("{0}\t", i);
                
            }

            file = @"ExamineZ.txt";
            using (StreamWriter output = File.AppendText(PolarBearal_OUTPUT_DIR + file))
            {
                for (int i = 0; i < examineZP.Length; i++)
                    output.Write("\n{0}\t{1}\t{2}\t{3}", _ZP[i], _ZN[i], _examineZP[i], _examineZN[i]);
            }



        }


        //Constructor 
        //uses a previously extracted barrel data to create a quick model of a barrel
        public PolarBearal(string PDBid)
        {
            // soluble barrels using 'all' method do not look at Z, but mono and poly do
            if (Global.METHOD == "all") { use_zone = false; }
            else { use_zone = true; }
            try
            {
                largestAltSeq = 0;
                proteinID = PDBid;
                totalProteins++;

                string file = @"betaBarrel_aaOnly/" + PDBid + ".txt";
                string[] FindChains = File.ReadAllLines(PolarBearal_OUTPUT_DIR + file);
                numStrands = FindChains.Length;

                simpleBarrel = new aa[numStrands][];
                Queue<aa> chain = new Queue<aa>();
                file = @"betaBarrelRawData/" + PDBid + ".txt";
                string[] BarrelChains = File.ReadAllLines(PolarBearal_OUTPUT_DIR + file);

                int curStrand = 0;
                char[] delPunc = { ':', ',' };
                foreach (string amino in BarrelChains)
                {
                    string[] preaa = amino.Split(delPunc);

                    if (preaa[0] != curStrand.ToString())
                    {
                        simpleBarrel[curStrand] = chain.ToArray();
                        chain.Clear();
                        curStrand++;
                    }

                    aa tempaa = new aa(preaa[1]);

                    if (preaa[2].ToLower() == "true")
                    {
                        tempaa.Inward = true;
                    }
                    else
                    {
                        tempaa.Inward = false;
                    }

                    tempaa.ResNum = Convert.ToInt16(preaa[3]);
                    tempaa.SeqID = Convert.ToInt16(preaa[4]);
                    tempaa.height = Double.Parse(preaa[5]);

                    tempaa.angle = Double.Parse(preaa[6]);
                    
                    chain.Enqueue(tempaa);
                }

                simpleBarrel[curStrand] = chain.ToArray();
                seperatedBarrel = new Queue<aa[][]>();
                SeperateBarrel();
                

                if (use_zone) { printZ(); }
                PolarBearalCalculations();
                PinNoutByProtein();
                //update all db counts since attempt has successfullly completed
                _totalProteins = totalProteins;
                _totalStrands = totalStrands;
                _totalAAs = totalAAs;
                _usedAAs = usedAAs;
                _seqCount = seqCount;
                _strandSizes = strandSizes;
                _AlternatingStrandSeqLengths = AlternatingStrandSeqLengths;
                _AlternatingInOutSeqLengths = AlternatingInOutSeqLengths;
                _AlternatingPolNonPSeqLengths = AlternatingPolNonPSeqLengths;
                _AlternatingBothSeqLengths = AlternatingBothSeqLengths;
                __aa = _aa;
                _aaInward = aaInward;
                _aaDevInward = aaDevInward;
                _aaOutward = aaOutward;
                _aaDevOutward = aaDevOutward;
                _aap1 = aap1;
                _aap2 = aap2;
                _aap3 = aap3;
                _aan1 = aan1;
                _aan2 = aan2;
                _aan3 = aan3;
                _ZP = ZP;
                _ZN = ZN;
                _examineZP = examineZP;
                _examineZN = examineZN;

            }
            catch
            {
                //revert all db counts to what they were at the end of last successful attempt
                totalProteins = _totalProteins;
                totalStrands = _totalStrands;
                totalAAs = _totalAAs;
                usedAAs = _usedAAs;
                seqCount = _seqCount;
                strandSizes = _strandSizes;
                AlternatingStrandSeqLengths = _AlternatingStrandSeqLengths;
                AlternatingInOutSeqLengths = _AlternatingInOutSeqLengths;
                AlternatingPolNonPSeqLengths = _AlternatingPolNonPSeqLengths;
                AlternatingBothSeqLengths = _AlternatingBothSeqLengths;
                _aa = __aa;
                aaInward = _aaInward;
                aaDevInward = _aaDevInward;
                aaOutward = _aaOutward;
                aaDevOutward = _aaDevOutward;
                aap1 = _aap1;
                aap2 = _aap2;
                aap3 = _aap3;
                aan1 = _aan1;
                aan2 = _aan2;
                aan3 = _aan3;
                ZP = _ZP;
                ZN = _ZN;
                examineZP = _examineZP;
                examineZN = _examineZN;
            }

        }

        //Constructor 
        //uses xml and the rest of the program to discover a barrel
        //extracts important data and add it to file for quicker extraction in future program runs
        public PolarBearal(ref Barrel myBarrel)
        {
            // soluble barrels using 'all' method do not look at Z, but mono and poly do
            if (Global.METHOD == "all") { use_zone = false; }
            else { use_zone = true; }



            proteinID = myBarrel.PdbName;
            largestAltSeq = 0;
            totalProteins++;



            string use_chain = myBarrel.Strands[0].ChainName;
            numStrands = 0;
            for (int curStrand = 0; curStrand < myBarrel.Strands.Count; curStrand++)
            {
                // only save data for one chain to avoid overcounting poly-barrels (except for 2 poly exceptions)
                if (myBarrel.PdbName.ToUpper() == "4TW1" & (myBarrel.Strands[curStrand].ChainName == "E" | myBarrel.Strands[curStrand].ChainName == "F") ) numStrands++;
                else if(myBarrel.PdbName.ToUpper() == "3B07" & (myBarrel.Strands[curStrand].ChainName == "A" | myBarrel.Strands[curStrand].ChainName == "B")) numStrands++;
                else if (use_chain == myBarrel.Strands[curStrand].ChainName) numStrands++;
            }
            //2D array to store data;[strand number] [amino acid number]
            simpleBarrel = new aa[numStrands][];

            //create simple barrel
            int strand_ctr = 0;
            foreach (Strand s in myBarrel)
            {
                // only save data for one chain to avoid overcounting poly-barrels
                if ((myBarrel.PdbName.ToUpper() == "4TW1" & (s.ChainName == "E" | s.ChainName == "F")) |
                    (myBarrel.PdbName.ToUpper() == "3B07" & (s.ChainName == "A" | s.ChainName == "B")) |
                    (use_chain == s.ChainName))
                {
                    simpleBarrel[strand_ctr] = new aa[s.Residues.Count()];
                    //add aa
                    int aa_cnt = 0;
                    foreach (Res cur_aa in s)
                    {
                        //add aa member variables here
                        simpleBarrel[strand_ctr][aa_cnt] = new aa(cur_aa.OneLetCode);
                        simpleBarrel[strand_ctr][aa_cnt].Inward = cur_aa.Inward;
                        simpleBarrel[strand_ctr][aa_cnt].ResNum = cur_aa.ResNum;
                        simpleBarrel[strand_ctr][aa_cnt].SeqID = cur_aa.SeqID;
                        simpleBarrel[strand_ctr][aa_cnt].X = cur_aa.BackboneCoords["CA"].X;
                        simpleBarrel[strand_ctr][aa_cnt].Y = cur_aa.BackboneCoords["CA"].Y;
                        simpleBarrel[strand_ctr][aa_cnt].height = cur_aa.BackboneCoords["CA"].Z;
                        simpleBarrel[strand_ctr][aa_cnt].chain = cur_aa.ChainName;

                        simpleBarrel[strand_ctr][aa_cnt].angle = SharedFunctions.AngleBetween(cur_aa.BackboneCoords["CA"] - ((cur_aa.BackboneCoords["N"] + cur_aa.BackboneCoords["C"]) / 2), myBarrel.Axis);
                        aa_cnt++;
                    }
                    strand_ctr++;
                }
            }

            string file = @"betaBarrel_aaOnly/" + proteinID + ".txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + file))
            {
                for (int curStrand = 0; curStrand < simpleBarrel.GetLength(0); curStrand++)

                {
                    //add aa
                    for (int cur_aa = 0; cur_aa < simpleBarrel[curStrand].Length; cur_aa++)
                    {
                        output.Write(" {0} ", simpleBarrel[curStrand][cur_aa].m_aa_ID);
                    }
                    output.WriteLine();
                }
            }

            file = @"betaBarrelRawData/" + proteinID + ".txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + file))
            {
                for (int curStrand = 0; curStrand < simpleBarrel.GetLength(0); curStrand++)
                {
                    //add aa
                    for (int cur_aa = 0; cur_aa < simpleBarrel[curStrand].Length; cur_aa++)
                    {
                        output.WriteLine("{0}:{1},{2},{3},{4},{5},{6}", curStrand, simpleBarrel[curStrand][cur_aa].m_aa_ID, simpleBarrel[curStrand][cur_aa].Inward, simpleBarrel[curStrand][cur_aa].ResNum, simpleBarrel[curStrand][cur_aa].SeqID, simpleBarrel[curStrand][cur_aa].height, simpleBarrel[curStrand][cur_aa].angle);
                    }
                }
            }


            string _pdb;
            int _pdb_strands;
            string _res;
            int _res_num, _seqID;
            int _res_strand;
            double _res_ca_x;
            double _res_ca_y;
            double _res_ca_z;
            bool _inward;
            double _angle;
            string _chain;
            file = @"betaBarrelStrands/" + proteinID + ".txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + file))
            {
                output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", "_pdb", "_pdb_strands", "_res", "_res_num", "_seqID", "_res_strand", "_res_ca_x", "_res_ca_y", "_res_ca_z", "_inward", "_angle", "_chain");
                for (int curStrand = 0; curStrand < simpleBarrel.GetLength(0); curStrand++)
                {
                    //add aa
                    for (int cur_aa = 0; cur_aa < simpleBarrel[curStrand].Length; cur_aa++)
                    {
                        _pdb = proteinID;
                        _pdb_strands = simpleBarrel.GetLength(0);
                        _res = simpleBarrel[curStrand][cur_aa].m_aa_ID;
                        _res_num = simpleBarrel[curStrand][cur_aa].ResNum + 1;
                        _seqID = simpleBarrel[curStrand][cur_aa].SeqID;
                        _res_strand = curStrand + 1;
                        _res_ca_x = simpleBarrel[curStrand][cur_aa].X;
                        _res_ca_y = simpleBarrel[curStrand][cur_aa].Y;
                        _res_ca_z = simpleBarrel[curStrand][cur_aa].height;
                        _inward = simpleBarrel[curStrand][cur_aa].Inward;
                        _angle = simpleBarrel[curStrand][cur_aa].angle;
                        _chain = simpleBarrel[curStrand][cur_aa].chain;
                        output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", _pdb, _pdb_strands, _res, _res_num, _seqID, _res_strand, _res_ca_x, _res_ca_y, _res_ca_z, _inward, _angle, _chain);
                    }
                }
            }
            // The following file can be used to only output residues that we want to use for calculations for soluble vs. membrane barrel paper
            // soluble = all strand residues
            // mono = strand residues within membrane region
            // poly = non-repeating strand residues within membrane region
            file = @"betaBarrelStrands_UseForCalc/" + proteinID + "_LimitedBarrelStrandResidues.txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(PolarBearal_OUTPUT_DIR + file))
            {
                output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", "_pdb", "_pdb_strands", "_res", "_res_num", "_seqID", "_res_strand", "_res_ca_x", "_res_ca_y", "_res_ca_z", "_inward", "_angle", "_chain");
                for (int curStrand = 0; curStrand < simpleBarrel.GetLength(0); curStrand++)
                {
                    _chain = simpleBarrel[curStrand][0].chain;
                    if ((myBarrel.PdbName.ToUpper() == "4TW1" & (_chain == "E" | _chain == "F")) |
                    (myBarrel.PdbName.ToUpper() == "3B07" & (_chain == "A" | _chain == "B")) |
                    (use_chain == _chain))
                    {
                        //add aa
                        for (int cur_aa = 0; cur_aa < simpleBarrel[curStrand].Length; cur_aa++)
                        {
                            _res_ca_z = simpleBarrel[curStrand][cur_aa].height;
                            if (!use_zone || (_res_ca_z > -zone && _res_ca_z < zone))
                            {
                                _pdb = proteinID;
                                _pdb_strands = simpleBarrel.GetLength(0);
                                _res = simpleBarrel[curStrand][cur_aa].m_aa_ID;
                                _res_num = simpleBarrel[curStrand][cur_aa].ResNum + 1;
                                _seqID = simpleBarrel[curStrand][cur_aa].SeqID;
                                _res_strand = curStrand + 1;
                                _res_ca_x = simpleBarrel[curStrand][cur_aa].X;
                                _res_ca_y = simpleBarrel[curStrand][cur_aa].Y;
                                _inward = simpleBarrel[curStrand][cur_aa].Inward;
                                _angle = simpleBarrel[curStrand][cur_aa].angle;
                                _chain = simpleBarrel[curStrand][cur_aa].chain;
                                output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", _pdb, _pdb_strands, _res, _res_num, _seqID, _res_strand, _res_ca_x, _res_ca_y, _res_ca_z, _inward, _angle, _chain);
                            }
                        }
                    }
                }
            }
            
        }


        //breaks a strand into substrands
        public void SeperateBarrel()
        {
            Queue<aa[]> seperatedChain = new Queue<aa[]>();
            int subChainCounter = 0;
            Queue<aa> subChain = new Queue<aa>();
            aa[] chain;

            for (int curStrand = 0; curStrand < numStrands; curStrand++)
            {
                chain = simpleBarrel[curStrand];

                //break apart chain
                for (int cur_aa = 0; cur_aa < chain.Length; cur_aa++)
                {
                    if (subChain.Count == 0)
                    {
                        subChain.Enqueue(chain[cur_aa]);
                    }
                    else if (!alternating_aa(chain[cur_aa - 1], chain[cur_aa]))
                    {
                        //check for an increased largest chain count
                        if (subChain.Count > largestAltSeq) largestAltSeq = subChain.Count;

                        aa[] subChain_aa = new aa[subChain.Count];
                        subChainCounter = 0;
                        foreach (aa number in subChain)
                        {
                            subChain_aa[subChainCounter++] = number;
                        }

                        subChain.Clear();
                        subChain.Enqueue(chain[cur_aa]);

                        seperatedChain.Enqueue(subChain_aa);
                    }
                    else
                    {
                        subChain.Enqueue(chain[cur_aa]);
                    }
                }
                //adds last chain if it ended with alternating sequances
                if (subChain.Count > 0)
                {
                    if (subChain.Count > largestAltSeq) largestAltSeq = subChain.Count;

                    aa[] subChain_aa = new aa[subChain.Count];
                    subChainCounter = 0;
                    foreach (aa number in subChain)
                    {
                        subChain_aa[subChainCounter++] = number;
                    }
                    subChain.Clear();
                    seperatedChain.Enqueue(subChain_aa);
                }


                aa[][] newChain = new aa[seperatedChain.Count][];
                for (int curSubChain = 0; 0 < seperatedChain.Count; curSubChain++)
                {
                    newChain[curSubChain] = seperatedChain.Dequeue();
                }

                seperatedBarrel.Enqueue(newChain);
            }
        }

        //determines if two given AA are PtoN or NtoP
        private bool alternating_aa(aa prevAA, aa curAA)
        {
            if (prevAA != null && curAA != null)
            {
                if ('P' == prevAA.mPolarity || prevAA.mPolarity == 'Q')
                {
                    if ('N' == curAA.mPolarity || curAA.mPolarity == 'Q') return (true);
                    else return (false);
                }
                else if ('N' == prevAA.mPolarity || prevAA.mPolarity == 'Q')
                {
                    if ('P' == curAA.mPolarity || curAA.mPolarity == 'Q') return (true);
                    else return (false);
                }
                else//no other way to be alternating
                {
                    return (false);
                }
            }
            return (false);
        }


        //Displays IN, Pin, OUT, Nout for each protein
        public void PinNoutByProtein()
        {
            int AAs = 0, Pin = 0, Nout = 0, IN = 0, OUT = 0;

            foreach (aa[] strand in this.simpleBarrel)
            {
                foreach (aa amino in strand)
                {
                    AAs++;

                    if (amino.Inward)
                    {
                        IN++;
                        if (amino.mPolarity == 'P' || amino.mPolarity == 'Q') Pin++;
                    }

                    if (!amino.Inward)
                    {
                        OUT++;
                        if (amino.mPolarity == 'N' || amino.mPolarity == 'Q') Nout++;
                    }
                }
            }

            string file = @"PinNoutByProtein.txt";
            using (System.IO.StreamWriter output = File.AppendText(PolarBearal_OUTPUT_DIR + file))
            {
                output.Write("\n {0} \t {1} \t {2} \t {3} \t {4} \t {5}", proteinID, AAs, IN, Pin, OUT, Nout);
            }
        }

        //For normal calculations
        public void PolarBearalCalculations()
        {
            totalStrands += simpleBarrel.Count();
            for (int i = 0; i < simpleBarrel.Length; i++)
            {
                strandSizes[simpleBarrel[i].Length]++;
            }

            foreach (aa[] strand in simpleBarrel)
            {
                foreach (aa amino in strand)
                {

                    if (!use_zone || (amino.height > -zone && amino.height < zone))
                    {
                        _aa[IntOfeRes(amino.m_aa_ID)]++;
                        if (amino.Inward) aaInward[IntOfeRes(amino.m_aa_ID)]++;
                        else aaOutward[IntOfeRes(amino.m_aa_ID)]++;
                    }
                }
            }

            //numbers to make a graph of occurances by chain length
            int[] OcurrancesPerLengthOfAltSeq = new int[largestAltSeq + 1];

            aa[][] temp2D_aa;
            int tempLength;
            for (int curStrand = 0; curStrand < numStrands; curStrand++)
            {
                temp2D_aa = seperatedBarrel.Dequeue();
                for (int curSubChain = 0; curSubChain < temp2D_aa.Length; curSubChain++)
                {
                    tempLength = temp2D_aa[curSubChain].Length;
                    OcurrancesPerLengthOfAltSeq[tempLength]++;
                    seqCount[tempLength]++;
                }
                seperatedBarrel.Enqueue(temp2D_aa);
            }

            string file = @"PolarBearalGraphOccurances.txt";
            using (StreamWriter sw = File.AppendText(PolarBearal_OUTPUT_DIR + file))
            {
                sw.Write("{0} \t", proteinID);
                foreach (int occurances in OcurrancesPerLengthOfAltSeq)
                {
                    sw.Write("\t {0}", occurances);
                }
                sw.Write("\n");
            }

            //AlternatingFacing();
            AlternatingSequenceLengths();
        }

        //converts amino acid 1 letter ID into associated interger
        public int IntOfeRes(string aa)
        {
            switch (aa)
            {
                case "L":
                    return (1);
                case "W":
                    return (2);
                case "P":
                    return (3);
                case "V":
                    return (4);
                case "I":
                    return (5);
                case "F":
                    return (6);
                case "Y":
                    return (7);
                case "H":
                    return (8);
                case "A":
                    return (9);
                case "M":
                    return (10);
                case "G":
                    return (11);
                case "T":
                    return (12);
                case "Q":
                    return (13);
                case "D":
                    return (14);
                case "N":
                    return (15);
                case "K":
                    return (16);
                case "S":
                    return (17);
                case "R":
                    return (18);
                case "E":
                    return (19);
                case "C":
                    return (20);
                default:
                    return (0);
            }
        }

        //use to examine specific values at different Zs
        public void printZ()
        {
            foreach (aa[] strand in simpleBarrel)
            {
                bool first = true;
                aa previousAA = strand[0];
                int strandNum = 0;
                foreach (aa amino in strand)
                {
                    if (!first)
                    {
                        if (amino.height > -zone && amino.height < zone && previousAA.height > -zone && previousAA.height < zone)
                        {
                            if (amino.mPolarity == 'P' || amino.mPolarity == 'Q') ZP[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            if (amino.mPolarity == 'N' || amino.mPolarity == 'Q') ZN[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;

                            if (amino.mPolarity == 'P' && previousAA.mPolarity == 'P')// || previousAA.mPolarity == 'Q' || amino.mPolarity == 'Q')
                            {
                                examineZP[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            }
                            else if (amino.mPolarity == 'N' && previousAA.mPolarity == 'N')// || previousAA.mPolarity == 'Q' || amino.mPolarity == 'Q')
                            {
                                examineZN[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            }
                            else { }


                        }
                        strandNum++;
                    }
                    first = false;
                    previousAA = amino;
                }
            }
        }

        //Examine Alternating Directionality
        void AlternatingFacing()
        {
            bool prevAAInward = true;
            double prevAAHeight = 0.0;
            char prevAAPol = 'P';
            int aaCount, seqCount, seqCount2;
            aa aa1_PreDev = new aa("A"), aa2_Dev = new aa("A"), aa3_PostDev = new aa("A");
            foreach (aa[] strand in this.simpleBarrel)
            {
                aaCount = 0;
                seqCount = 0;
                seqCount2 = 0;
                bool range = false;
                foreach (aa amino in strand)
                {
                    if (!use_zone || (amino.height > -zone && amino.height < zone) )
                    {
                        totalAAs++;
                        range = true;
                        if (aaCount > 0)
                        {
                            //totalAAs++;
                            if (seqCount == 0)
                            {
                                seqCount++;
                            }
                            else
                            {
                                if (amino.Inward != prevAAInward) seqCount++;
                                else
                                {
                                    AlternatingInOutSeqLengths[seqCount]++;
                                    if (seqCount != 1)
                                    {
                                        if (amino.Inward && prevAAInward) aaDevInward[IntOfeRes(amino.m_aa_ID)]++;
                                        else if (!amino.Inward && !prevAAInward) aaDevOutward[IntOfeRes(amino.m_aa_ID)]++;
                                    }
                                    seqCount = 1;
                                }



                                if (seqCount2 == 0)
                                {
                                    seqCount2++;
                                    aa3_PostDev = amino;
                                }
                                else if (seqCount2 == 1)
                                {
                                    seqCount2++;
                                    aa2_Dev = aa3_PostDev;
                                    aa3_PostDev = amino;
                                }
                                else
                                {
                                    seqCount2++;
                                    aa1_PreDev = aa2_Dev;
                                    aa2_Dev = aa3_PostDev;
                                    aa3_PostDev = amino;

                                    if (aa1_PreDev.mPolarity == 'P' && aa2_Dev.mPolarity == 'P')
                                    {
                                        aap1[IntOfeRes(aa1_PreDev.m_aa_ID)]++;
                                        aap2[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                                        aan3[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                                    }

                                    if (aa1_PreDev.mPolarity == 'N' && aa2_Dev.mPolarity == 'N')
                                    {
                                        aan1[IntOfeRes(aa1_PreDev.m_aa_ID)]++;
                                        aan2[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                                        aap3[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                                    }
                                }
                            }
                        }
                        else if (use_zone){ seqCount++;}
                    }


                    prevAAInward = amino.Inward;
                    prevAAHeight = amino.height;
                    prevAAPol = amino.mPolarity;
                    aaCount++;
                }
                if (aa3_PostDev.mPolarity == 'P' && aa2_Dev.mPolarity == 'P')
                {
                    aap1[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                    aap2[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                    aan3[21]++;
                }

                if (aa3_PostDev.mPolarity == 'N' && aa2_Dev.mPolarity == 'N')
                {
                    aan1[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                    aan2[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                    aap3[21]++;
                }
                if (range) AlternatingInOutSeqLengths[seqCount]++;
            }
        }

        // get lengths of alternating direction and hydrophobicity seq. (new cleaner version of AlternatingFacing
        void AlternatingSequenceLengths()
        {
            bool prevAAInward = true;
            double prevAAHeight = 0.0;
            char prevAAPol = 'P';
            int numAllStrandAAs, numStrandAAs, altDirSeqLen, altHydSeqLen, altBothSeqLen;

            int strand_ctr = 0;
            foreach (aa[] strand in this.simpleBarrel)
            {
                numStrandAAs = 0; numAllStrandAAs = 0;
                altDirSeqLen = 0;
                altHydSeqLen = 0;
                altBothSeqLen = 0;
                
                foreach (aa amino in strand)
                {
                    totalAAs++; numAllStrandAAs++;
                    if (!use_zone || (amino.height > -zone && amino.height < zone))
                    {
                        usedAAs++;
                        if (numStrandAAs == 0)
                        {
                            altDirSeqLen++;
                            altHydSeqLen++;
                            altBothSeqLen++;
                        }
                        else
                        {
                            // count alternating direction sequences
                            if (amino.Inward != prevAAInward) { altDirSeqLen++; }
                            else
                            {
                                // failed to alternate, so add previous alternating direction sequence to static counts and reset 
                                AlternatingInOutSeqLengths[altDirSeqLen]++;
                                altDirSeqLen = 1;
                            }

                            // count alternating hydrophobicity sequences
                            if (amino.mPolarity != prevAAPol) { altHydSeqLen++; }
                            else
                            {
                                // failed to alternate, so add previous alternating direction sequence to static counts and reset 
                                AlternatingPolNonPSeqLengths[altHydSeqLen]++;
                                altHydSeqLen = 1;
                            }

                            // count alternating direction hydrophobicity sequences
                            if ((amino.mPolarity != prevAAPol) && (amino.Inward != prevAAInward)) { altBothSeqLen++; }
                            else
                            {
                                // failed to alternate, so add previous alternating direction sequence to static counts and reset 
                                AlternatingBothSeqLengths[altBothSeqLen]++;
                                altBothSeqLen = 1;
                            }
                        }

                    

                    prevAAInward = amino.Inward;
                    prevAAHeight = amino.height;
                    prevAAPol = amino.mPolarity;
                    numStrandAAs++;
                }
                }
                AlternatingInOutSeqLengths[altDirSeqLen]++;
                AlternatingPolNonPSeqLengths[altHydSeqLen]++;
                AlternatingBothSeqLengths[altBothSeqLen]++;
                AlternatingStrandSeqLengths[numStrandAAs]++;
                if (0 == numStrandAAs)
                {
                    Console.WriteLine("PDB {0} strand {1} tried to add empty sequence. Strand has {2} residues, {3} are in zone!!!", this.proteinID, strand_ctr, numAllStrandAAs, numStrandAAs);
                }
                strand_ctr++;
            }
        }

        //member class aa
        class aa
        {
            //contructors
            public aa(string aa_ID)
            {
                m_aa_ID = aa_ID;
                mPolarity = 'U';
                assignPolarity();
                Inward = false;
                ResNum = -100;
                SeqID = -100;
                height = -100;
                chain = "";
            }

            //assigns polarity from amino acid ID
            public void assignPolarity()
            {
                string ID = this.m_aa_ID;

                if (P.Contains(ID))
                {
                    mPolarity = 'P';
                }
                else if (N.Contains(ID))
                {
                    mPolarity = 'N';
                }
                else
                {
                    mPolarity = 'Q';
                }
            }

            //member variables
            public string m_aa_ID;//stores amino acid's 3 letter abreviation
            public string chain;// stores chainID for this aa
            public char mPolarity;//should be P,N, or U(represent unassigned polarity}
            public bool Inward;//is true if amino acid is facing inward, false if facing outward
            public int ResNum;//residue number in protein
            public int SeqID;//residue sequence position (includes unresolved residues and may not start at 1)
            public double X;//height along z axis
            public double Y;//height along z axis
            public double height;//height along z axis
            public double angle;
        }


        //member variables
        string proteinID;
        aa[][] simpleBarrel;
        Queue<aa[][]> seperatedBarrel;
        int numStrands, largestAltSeq;
        // probably a better way to do this, but i duplicated all the following variables so that they
        // can be reverted to previous successful attempt should current attempt fail (so only working barrels in db are counted)
        public static int totalProteins = 0, totalStrands = 0, totalAAs = 0, usedAAs=0;
        public static int _totalProteins = 0, _totalStrands = 0, _totalAAs = 0, _usedAAs=0;

        public static int[] seqCount = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        public static int[] _seqCount = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        public static string P = "DEKNQRSTCH", N = "AFILVWYMP", Q = "G";

        public static int[] _aa = new int[21];
        public static int[] __aa = new int[21];

        public static int[] aaInward = new int[21];
        public static int[] _aaInward = new int[21];

        public static int[] aaOutward = new int[21];
        public static int[] _aaOutward = new int[21];

        public static int[] aaDevInward = new int[21];
        public static int[] _aaDevInward = new int[21];

        public static int[] aaDevOutward = new int[21];
        public static int[] _aaDevOutward = new int[21];

        public static int[] strandSizes = new int[50];
        public static int[] _strandSizes = new int[50];
        public static int[] ZP = new int[150];
        public static int[] _ZP = new int[150];
        public static int[] ZN = new int[150];
        public static int[] _ZN = new int[150];
        public static int[] AlternatingInOutSeqLengths = new int[100];
        public static int[] _AlternatingInOutSeqLengths = new int[100];
        public static int[] AlternatingPolNonPSeqLengths = new int[100];
        public static int[] _AlternatingPolNonPSeqLengths = new int[100];
        public static int[] AlternatingBothSeqLengths = new int[100];
        public static int[] _AlternatingBothSeqLengths = new int[100];
        public static int[] AlternatingStrandSeqLengths = new int[100];
        public static int[] _AlternatingStrandSeqLengths = new int[100];
        public static int[] examineZP = new int[125];
        public static int[] _examineZP = new int[125];
        public static int[] examineZN = new int[125];
        public static int[] _examineZN = new int[125];
        public static int numDoubleChecks = 0;

        //new data 10-3-15 Ryan
        public static int[] aap1 = new int[21];
        public static int[] aan1 = new int[21];
        public static int[] aap2 = new int[21];
        public static int[] aan2 = new int[21];
        public static int[] aap3 = new int[22];
        public static int[] aan3 = new int[22];

        public static int[] _aap1 = new int[21];
        public static int[] _aan1 = new int[21];
        public static int[] _aap2 = new int[21];
        public static int[] _aan2 = new int[21];
        public static int[] _aap3 = new int[22];
        public static int[] _aan3 = new int[22];

        public static bool use_zone = true;
    }
}