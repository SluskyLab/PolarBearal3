/*
**  File: Program.cs
**  Started: pre 6/01/2015
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
using betaBarrelProgram.AtomParser;
using betaBarrelProgram.BarrelStructures;
using betaBarrelProgram.Mono;


using System.Collections;
using System.IO;

using System.Net.Http;// for downloading new PDBs 


namespace betaBarrelProgram
{
    public static class Global
    {
        public static string POLARBEARAL_DIR = Directory.GetParent(System.IO.Directory.GetCurrentDirectory()).Parent.Parent.Parent.FullName;

        public static string MONO_DB_DIR = POLARBEARAL_DIR + "./DB/PDBs/";
        public static string MONO_DB_file = POLARBEARAL_DIR + "./DB/MonoDBList.txt";
        public static string MONO_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/mono/";

        public static string POLY_DB_DIR = POLARBEARAL_DIR + "./DB/PolyBarrelsDB/";
        public static string POLY_DB_file = POLARBEARAL_DIR + "./DB/PolyDBList.txt";
        public static string POLY_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/poly/";

        public static string MEMB_DB_DIR = POLARBEARAL_DIR + "./DB/PDBs/";
        public static string MEMB_DB_file = POLARBEARAL_DIR + "./DB/MembDBList.txt";
        public static string MEMB_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/membrane/";

        public static string SOLUBLE_DB_DIR = POLARBEARAL_DIR + "./DB/PDBs/";
        public static string SOLUBLE_DB_file = POLARBEARAL_DIR + "./DB/SolubleDBList.txt";
        public static string SOLUBLE_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/soluble/";

        // change method and TEST_DB_DIR to test subsets of PDBs with mono, poly, and all methods
        public static string TEST_DB_file = POLARBEARAL_DIR + "./DB/TestDBList.txt";
        public static string TEST_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/test/";
        public static string TEST_DB_DIR = POLARBEARAL_DIR + "./DB/AF_PDBs/";
        public static string TEST_METHOD = "mono";
        // default to comprehensive input and output
        public static string DB_DIR = POLARBEARAL_DIR + "./DB/PDBs/";
        public static string OUTPUT_DIR = POLARBEARAL_DIR + "./Output/";
        public static string DB_file = POLARBEARAL_DIR + "./DB/AllDBList.txt";
        public static string METHOD = POLARBEARAL_DIR + "mono";

        public static string AF_DB_DIR = POLARBEARAL_DIR + "./DB/AF_PDBs/";
        public static string AF_DB_file = POLARBEARAL_DIR + "./DB/AF_DBList.txt";
        public static string AF_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/AF_mono/";

        public static string AF2_DB_file = POLARBEARAL_DIR + "./DB/AF_DBList2.txt";
        public static string AF2_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/AF_mono_test/";

        public static string AFWeird_DB_file = POLARBEARAL_DIR + "./DB/AF_DBweird.txt";
        public static string AFWeird_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/AF_weird_test/";

        public static string AF3_DB_file = POLARBEARAL_DIR + "./DB/AF_DBList3.txt";
        public static string AF3_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/AF_mono_test2/";

        public static string BIG1_DB_DIR = "./DB/big_af/big_af_set_PDBs/top_preds/";
        public static string BIG1_DB_file = "./DB/big_af/big_af_set_PDBs/top_preds/quality80_pdb_list.txt";
        public static string BIG1_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/big_test1/";

        public static string CUSTOM_DB_DIR = "./DB/CustomPDBs/";
        public static string CUSTOM_DB_file = "./DB/CustomDBList.txt";
        public static string CUSTOM_OUTPUT_DIR = POLARBEARAL_DIR + "./Output/custom/";

        public static void change_to_BIG_data()
        {
            DB_DIR = BIG1_DB_DIR;
            DB_file = BIG1_DB_file;
            OUTPUT_DIR = BIG1_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_AF3_data()
        {
            DB_DIR = AF_DB_DIR;
            DB_file = AF3_DB_file;
            OUTPUT_DIR = AF3_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_AFWeird_data()
        {
            DB_DIR = AF_DB_DIR;
            DB_file = AFWeird_DB_file;
            OUTPUT_DIR = AFWeird_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_AF2_data()
        {
            DB_DIR = AF_DB_DIR;
            DB_file = AF2_DB_file;
            OUTPUT_DIR = AF2_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_AF_data()
        {
            DB_DIR = AF_DB_DIR;
            DB_file = AF_DB_file;
            OUTPUT_DIR = AF_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_mono_data()
        {
            DB_DIR = MONO_DB_DIR;
            DB_file = MONO_DB_file;
            OUTPUT_DIR = MONO_OUTPUT_DIR;
            METHOD = "mono";
        }
        public static void change_to_poly_data()
        {
            DB_DIR = POLY_DB_DIR;
            DB_file = POLY_DB_file;
            OUTPUT_DIR = POLY_OUTPUT_DIR;
            METHOD = "poly";
        }
        public static void change_to_membrane_data()
        {
            DB_DIR = MEMB_DB_DIR;
            DB_file = MEMB_DB_file;
            OUTPUT_DIR = MEMB_OUTPUT_DIR;
            METHOD = "mono/poly";
        }
        public static void change_to_soluble_data()
        {
            DB_DIR = SOLUBLE_DB_DIR;
            DB_file = SOLUBLE_DB_file;
            OUTPUT_DIR = SOLUBLE_OUTPUT_DIR;
            METHOD = "all";
        }

        public static void change_to_test_data()
        {
            DB_DIR = TEST_DB_DIR;
            DB_file = TEST_DB_file;
            OUTPUT_DIR = TEST_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_to_custom_data()
        {
            DB_DIR = CUSTOM_DB_DIR;
            DB_file = CUSTOM_DB_file;
            OUTPUT_DIR = CUSTOM_OUTPUT_DIR;
            METHOD = "mono";
        }

        public static void change_dataset()
        {
            Console.WriteLine("0. mono");
            Console.WriteLine("1. poly (obsolete)");
            Console.WriteLine("2. membrane");
            Console.WriteLine("3. soluble");
            Console.WriteLine("4. test (multi-method, in development)");
            Console.WriteLine("5. AF ()");
            Console.WriteLine("6. AF test (obsolete)");
            Console.WriteLine("7. AF weird (obsolete)");
            Console.WriteLine("8. AF test2 (obsolete)");
            Console.WriteLine("9. big testing 1 (obsolete)");
            Console.WriteLine("10. custom");

            string user_input = "";
            user_input = Console.ReadLine();

            switch (user_input)
            {
                case "0":
                    change_to_mono_data();
                    break;
                case "1":
                    change_to_poly_data();
                    break;
                case "2":
                    change_to_membrane_data();
                    break;
                case "3":
                    change_to_soluble_data();
                    break;
                case "4":
                    change_to_test_data();
                    break;
                case "5":
                    change_to_AF_data();
                    break;
                case "6":
                    change_to_AF2_data();
                    break;
                case "7":
                    change_to_AFWeird_data();
                    break;
                case "8":
                    change_to_AF3_data();
                    break;
                case "9":
                    change_to_BIG_data();
                    break;
                case "10":
                    change_to_custom_data();
                default:
                    change_to_mono_data();
                    break;
            }
            SharedFunctions.create_dir(OUTPUT_DIR);
        }


        public static string parameterFile = POLARBEARAL_DIR + "PolarBearal/par_hbond_1.txt";
        public static Dictionary<string, AminoAcid> AADict = SharedFunctions.makeAADict();

        //The values in this dictionary are transcribed from CHARMM36 all-hydrogen topology file for proteins, May 2011
        public static Dictionary<Tuple<string, string>, double> partialChargesDict = new Dictionary<Tuple<string, string>, double>
                {
                    {new Tuple<string, string>("ARG", "HE"), 0.44},
                    {new Tuple<string, string>("ARG", "NE"), -0.70},
                    {new Tuple<string, string>("ARG", "NH1"), -0.80},
                    {new Tuple<string, string>("ARG", "HH11"), 0.46},
                    {new Tuple<string, string>("ARG", "HH12"), 0.46},
                    {new Tuple<string, string>("ARG", "NH2"), -0.80},
                    {new Tuple<string, string>("ARG", "HH21"), 0.46},
                    {new Tuple<string, string>("ARG", "HH22"), 0.46},

                    {new Tuple<string, string>("ASN", "ND2"), -0.62},
                    {new Tuple<string, string>("ASN", "HD21"), 0.32},
                    {new Tuple<string, string>("ASN", "HD22"), 0.40},
                    {new Tuple<string, string>("ASN", "CG"), 0.55},
                    {new Tuple<string, string>("ASN", "OD1"), -0.55},

                    {new Tuple<string, string>("ASP", "CG"), 0.62},
                    {new Tuple<string, string>("ASP", "OD1"), -0.76},
                    {new Tuple<string, string>("ASP", "OD2"), -0.76},

                    {new Tuple<string, string>("CYS", "HG"), 0.16},
                    {new Tuple<string, string>("CYS", "SG"), -0.23},

                    {new Tuple<string, string>("GLN", "NE2"), -0.62},
                    {new Tuple<string, string>("GLN", "HE21"), 0.32},
                    {new Tuple<string, string>("GLN", "HE22"), 0.30},
                    {new Tuple<string, string>("GLN", "CD"), 0.55},
                    {new Tuple<string, string>("GLN", "OE1"), -0.55},

                    {new Tuple<string, string>("GLU", "CD"), 0.62},
                    {new Tuple<string, string>("GLU", "OE1"), -0.76},
                    {new Tuple<string, string>("GLU", "OE2"), -0.76},

                    {new Tuple<string, string>("GLY", "O"), -0.51},

                    {new Tuple<string, string>("HIS", "NE2"), -0.36},
                    {new Tuple<string, string>("HIS", "HE2"), 0.32},
                    {new Tuple<string, string>("HIS", "ND1"), -0.70},

                    {new Tuple<string, string>("LYS", "NZ"), -0.30},
                    {new Tuple<string, string>("LYS", "HZ1"), 0.33},
                    {new Tuple<string, string>("LYS", "HZ2"), 0.33},
                    {new Tuple<string, string>("LYS", "HZ3"), 0.33},

                    {new Tuple<string, string>("SER", "HG1"), 0.43},
                    {new Tuple<string, string>("SER", "OG1"), -0.66},

                    {new Tuple<string, string>("THR", "HG1"), 0.43},
                    {new Tuple<string, string>("THR", "OG1"), -0.66},

                    {new Tuple<string, string>("TRP", "NE1"), -0.51},
                    {new Tuple<string, string>("TRP", "HE1"), 0.37},

                    {new Tuple<string, string>("TYR", "HH"), 0.43},
                    {new Tuple<string, string>("TYR", "OH"), -0.54}

                };
    }

    class Program
    {
        static public void Display_menu()
        {
            Console.WriteLine("1. Test Structure for β-Barrel Validity"); 
            Console.WriteLine("2. Test Structure for Single-Barrel Analysis");
            Console.WriteLine("3. Choose dataset to run mass analysis (defaults to mono method)");
            Console.WriteLine("4. Run beta-barrel analysis for selected database");
            Console.WriteLine("5. PolarBearel"); //the polar bear
            Console.WriteLine("6. Multi-method analysis of a structure (in development)");
            Console.WriteLine("7. Run centroid analysis");
            Console.WriteLine("8. Run high-low analysis");
            Console.WriteLine("9. Print information about PDB files");
            Console.WriteLine("10. Run β-barrel pocket analysis for database structures");
            Console.WriteLine("Note: any other input will terminate this session. Please enter a number from 1-10.");

        }

        static void Main(string[] args)
        {
            DateTime startTime = DateTime.Now;
            //using (StreamWriter log = File.AppendText(Global.OUTPUT_DIR + "log.txt"))
            //{
            //    log.WriteLine("working with db ({0})", DateTime.Now);
            //    //log.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", "PdbName", "Success", "barrelType", "chainID", "numberOfStrands");
            //}

            string choice = "";

            string method_input = "3";
            string use_method = "all";

            while (choice != "10")
            {
                Display_menu();
                choice = Console.ReadLine();
                switch (choice)
                {
                    case "1":
                        // Test B-Barrel: check if a specific PDB can have its barrel read in and works
                        BarrelStructures.Protein _protein = null;
                        BarrelStructures.Barrel _barrel = null;
                        Console.WriteLine("Enter pdb:");
                        string PDBid = Console.ReadLine();

                        Console.WriteLine("Enter method (1=mono, 2=poly, default=all):");
                        method_input = Console.ReadLine();
                        use_method = "all";
                        if (method_input == "1"){ use_method = "mono"; }
                        if (method_input == "2") { use_method = "poly"; }
                        Console.WriteLine("Attempting to run {0} with {1} code.  Output will be in \n{2}", PDBid, use_method, Global.OUTPUT_DIR);

                        // change method from mono depending on wht method you want to use
                        runThisBetaBarrel(PDBid, use_method, ref _barrel, ref _protein);
                        break;
                    case "2":
                        BarrelEllipse.testEllipseSinglePDB();
                        //BarrelEllipse.runEllipseData();
                        break;
                    case "3":
                        Global.change_dataset();
                        Console.WriteLine("Using dataset in the file: {0}", Global.DB_file);
                        Console.WriteLine("Reading input files from: {0}", Global.DB_DIR);
                        Console.WriteLine("Writing output files to: {0}", Global.OUTPUT_DIR);
                        using (StreamWriter log = File.AppendText(Global.OUTPUT_DIR + "log.txt"))log.WriteLine("working with db ({0})", DateTime.Now);
                        break;
                    case "4":
                        RunBetaBarrelDatabase();
                        break;
                    case "5":
                        PolarBearal.RunPolarBearal();
                        break;
                    case "6":
                        Console.WriteLine("This has not been checked in a bit.  Proceed with caution");

                        Console.WriteLine("Enter method (1=mono, 2=poly, default=all):");
                        method_input = Console.ReadLine();
                        use_method = "all";
                        if (method_input == "1") { use_method = "mono"; }
                        if (method_input == "2") { use_method = "poly"; }
                        Console.WriteLine("Using {1} code.  Output will be in \n{2}", use_method, Global.OUTPUT_DIR);
                        SharedFunctions.RunCbeta2Axis(use_method);
                        break;
                    case "7":
                        Console.WriteLine("This has not been checked in a bit.  Proceed with caution");
                        SharedFunctions.run_Centroids();
                        break;
                    case "8":
                        Console.WriteLine("This has not been checked in a bit.  Proceed with caution");
                        BarrelEllipse.run_HighLowData();
                        break;
                    case "9":
                        PDBInfo testInfo = new PDBInfo();
                        break;
                    case "10":
                        SharedFunctions.run_DataForOMBBPockets();
                        break;
                    //case "":
                        //AlignPDB testAlign = new AlignPDB();
                    //    break;
                    default:
                        using (StreamWriter log = File.AppendText(Global.OUTPUT_DIR + "log.txt"))
                        {
                            log.WriteLine("ended program at {0} \n\n", DateTime.Now);
                        }
                        //10. Quit
                        choice = "10";
                        break;
                }
            }


            Console.WriteLine("Program took from: \n start: {0} \n   end: {1}", startTime, DateTime.Now);

            return;
        }


        public static Barrel runThisBetaBarrel(string pdb, string method)
        {
            //Protein ThisProtein = null;
            Barrel ThisBarrel = null;
            string PDB = pdb.Substring(0, 4).ToUpper();
            if (pdb.Count() == 6)
            {
                PDB = pdb.Substring(0, 6).ToUpper();
            }
            else if (pdb.Count() == 21)
            {
                PDB = pdb.Substring(0, 21);
            }
            else if (pdb.Count() == 25)
            {
                PDB = pdb.Substring(0, 25);
            }

            string pdbFileName = Global.DB_DIR + pdb + ".pdb";

            if (File.Exists(pdbFileName))
            {
                AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                //Console.WriteLine("opened {0}", pdbFileName);
                myAtomCat = Program.ReadPdbFile(pdbFileName, ref Global.partialChargesDict);
                Console.WriteLine("\nAttempting {0}", pdb);
                // use poly method
                if ("poly" == method)
                {
                    Console.WriteLine("Failed to run {0} using the poly methods\n", pdb);
                }
                // use all method 
                else if ("all" == method)
                {
                    Console.WriteLine("Failed to run {0} using the all methods\n", pdb);
                }
                // default to mono method
                else //if ("mono" == method)
                {
                    if(File.Exists(pdbFileName))
                    {
                        Console.WriteLine("Generating mono protein");
                        Protein newProt = new MonoProtein(ref myAtomCat, 0, PDB);
                        Console.WriteLine("Generating mono barrel");
                        ThisBarrel = new MonoBarrel(newProt.Chains[0], newProt);
                    }
                    else//catch
                    {
                        Console.WriteLine("Failed to run {0} using the mono methods\n", pdb);
                    }
                }
            }
            else
            {
                Console.WriteLine("could not find {0}", pdbFileName);
            }

            return (ThisBarrel);
        }

        public static void runThisBetaBarrel(string pdb, string method, ref Barrel ThisBarrel, ref Protein newProt)
        {
            string PDB = pdb;//.Substring(0, 4).ToUpper();
            string pdbFileName = Global.DB_DIR + pdb + ".pdb";

            if (File.Exists(pdbFileName))
            {
                AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                //Console.WriteLine("opened {0}", pdbFileName);
                myAtomCat = Program.ReadPdbFile(pdbFileName, ref Global.partialChargesDict);
                Console.WriteLine("\nAttempting {0}", pdb);
                // use poly method
                if ("poly" == method)
                {
                    Console.WriteLine("Failed to run {0} using the poly methods\n", pdb);
                }
                // use all method 
                else if ("all" == method)
                {
                   
                    Console.WriteLine("Failed to run {0} using the all methods\n", pdb);
                }
                // default to mono method
                else //if ("mono" == method)
                {
                    try
                    {
                        Console.WriteLine("Generating mono protein");
                        newProt = new MonoProtein(ref myAtomCat, 0, PDB);
                        Console.WriteLine("Generating mono barrel");
                        ThisBarrel = new MonoBarrel(newProt.Chains[0], newProt);
                    }
                    catch
                    {
                        Console.WriteLine("Failed to run {0} using the mono methods\n", pdb);
                    }
                }
            }
            else
            {
                Console.WriteLine("could not find {0}", pdbFileName);
            }
        }

        static public void RunBetaBarrelDatabase()
        {
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();

            string fileOfPDBs = Global.DB_file;
            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    string fileLocation2 = Global.OUTPUT_DIR + "AllBarrelChar.txt";
                    using (System.IO.StreamWriter AllBarrel_output = new System.IO.StreamWriter(fileLocation2))
                    {
                        //string newLine = "PDB" + "\t\t" + "Total Strands" +"\t" + "Length" + "\t" + "AvgLength" + "\t" + "MinLength" + "\t" + "MaxLength" + "\t" + "Radius" + "\t" + "Barrel Tilt";
                        //file.WriteLine(newLine);
                        AllBarrel_output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", "PDB", "total_strands", "strand_cnt", "avg_length", "min_length", "max_length", "avg_radius", "avg_tilt", "shear_num");
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

                                AllBarrel_output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", PDB, total_strands, strand_cnt, avg_length, min_length, max_length, avg_radius, avg_tilt, shear_num);


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

        //extracts atoms from pdb file
        public static AtomCategory ReadPdbFile(string pdbfilename, ref Dictionary<Tuple<string, string>, double> partialCharges)
        {
            AtomCategory myAtomCategory = new AtomCategory();

            ChainAtoms chainAtoms = new ChainAtoms();
            ArrayList atomList = new ArrayList();

            string preAsymId = "";
            string asymID = "";

            if (File.Exists(pdbfilename))
            {
                using (StreamReader sr = new StreamReader(pdbfilename))
                {
                    String line;
                    // Read and display lines from the file until the end of
                    // the file is reached.
                    Int32 seqID = -1;
                    Int32 previousSeqID = -1;
                    while ((line = sr.ReadLine()) != null)
                    {

                        string[] spLine = Array.FindAll<string>(((string)line).Split(
                        new char[] { ' ', '\t', ',', ';' }), delegate (string s) { return !String.IsNullOrEmpty(s); });
                        List<string> splitLine = new List<string>(spLine);
                        if (splitLine[0] == "HETATM" && splitLine[3] == "MSE")
                        {
                            splitLine[0] = "ATOM";
                            splitLine[3] = "MET";
                        }
                        if (splitLine[0] == "HETATM" && splitLine[3] == "CSO")
                        {
                            splitLine[0] = "ATOM";
                            splitLine[3] = "CYS";
                        }
                        if (splitLine[0] == "ENDMDL") break;

                        if (splitLine[0] == "ATOM")
                        {
                            AtomInfo myAtom = new AtomInfo();
                            asymID = line.Substring(21, 1);
                            myAtom.altConfID = "";
                            myAtom.atomId = Convert.ToInt32(line.Substring(6, 5).Trim());
                            myAtom.atomName = line.Substring(12, 4).Trim();
                            myAtom.residue = line.Substring(17, 3);
                            myAtom.authResidue = myAtom.residue;
                            myAtom.authSeqId = line.Substring(22, 4).Trim();
                            if (Convert.ToInt32(myAtom.authSeqId) != previousSeqID) seqID = seqID + 1;
                            myAtom.seqId = seqID.ToString();
                            previousSeqID = Convert.ToInt32(myAtom.authSeqId);

                            AtomParser.Coordinate coord = new Coordinate();
                            coord = new Coordinate(Convert.ToDouble(line.Substring(30, 8).Trim()), Convert.ToDouble(line.Substring(38, 8).Trim()), Convert.ToDouble(line.Substring(46, 8).Trim()));
                            myAtom.xyz = coord;
                            myAtom.atomType = line.Substring(76, 2).Trim();
                            //This handful of conditions are occasionally different and throw off the check for hydrogen bonding, etc
                            if (myAtom.atomName == "HG") myAtom.atomName = "HG1";
                            if (myAtom.atomName == "OG") myAtom.atomName = "OG1";
                            if ((myAtom.atomName == "NH1" && myAtom.residue == "ARG") || (myAtom.atomName == "NZ" && myAtom.residue == "LYS")) myAtom.atomType = "N1+";
                            if ((myAtom.atomName == "OD2" && myAtom.residue == "ASP") || (myAtom.atomName == "OE2" && myAtom.residue == "GLU")) myAtom.atomType = "O1-";

                            myAtom.bFac = Convert.ToDouble(line.Substring(60, 6).Trim());
                            myAtom.occupancy = 1.00;

                            //Compare chain ID. If this is now a new chain, write existing list, and create empty new ones.
                            if (preAsymId != asymID && preAsymId != "" && atomList.Count > 0)
                            {
                                chainAtoms.AsymChain = preAsymId;
                                chainAtoms.AuthAsymChain = preAsymId;
                                chainAtoms.EntityID = "1"; // problem with int to string in new version
                                chainAtoms.PolymerType = "-";

                                AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                                atomList.CopyTo(atomArray);
                                chainAtoms.CartnAtoms = atomArray;
                                myAtomCategory.AddChainAtoms(chainAtoms);
                                atomList = new ArrayList();
                                chainAtoms = new ChainAtoms();
                            }
                            //Add atom to current list of atoms if it is not a water molecule
                            if (myAtom.residue.ToUpper() != "HOH")
                            {
                                if (myAtom.atomType == "H")
                                {
                                    Tuple<string, string> key = new Tuple<string, string>(myAtom.residue, myAtom.atomName);
                                    if (partialCharges.ContainsKey(key) == true) atomList.Add(myAtom);
                                }
                                else atomList.Add(myAtom);
                            }
                            preAsymId = asymID;


                        }
                    }
                    //Capture final chain
                    if (atomList.Count > 0)
                    {
                        chainAtoms.AsymChain = preAsymId;
                        chainAtoms.AuthAsymChain = preAsymId;
                        chainAtoms.EntityID = "1";
                        chainAtoms.PolymerType = "-";
                        AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                        atomList.CopyTo(atomArray);
                        chainAtoms.CartnAtoms = atomArray;
                        myAtomCategory.AddChainAtoms(chainAtoms);
                        atomList = new ArrayList();
                        chainAtoms = new ChainAtoms();
                    }

                }
            }

            myAtomCategory.Resolution = 2.5;

            return myAtomCategory;
        }

        #region externalClassDefs

        public class TurnListSorter : IComparer<List<string>>
        {
            public int Compare(List<string> _l1, List<string> _l2)
            {
                if (_l1.Count != _l2.Count)
                { return _l2.Count.CompareTo(_l1.Count); } // returns in descending order
                return _l1[0].CompareTo(_l2[0]); // returns in alphabetical order if equal counts
            }
        }
        public class LKandCLIDSorter : IComparer<string> // format: loopKey_clusterID
        {
            public int Compare(string _s1, string _s2)
            {
                string lk1 = _s1.Substring(0, _s1.IndexOf("_"));
                int clID1 = Convert.ToInt32(_s1.Substring(_s1.IndexOf("_") + 1));
                string lk2 = _s2.Substring(0, _s2.IndexOf("_"));
                int clID2 = Convert.ToInt32(_s2.Substring(_s2.IndexOf("_") + 1));
                if (lk1 != lk2)
                { return lk1.CompareTo(lk2); }
                else
                { return clID1.CompareTo(clID2); }
            }
        }

        public class SupplDataSorter : IComparer<string>
        {
            public int Compare(string _s1, string _s2)
            {
                string[] splitIntoWords1 = Array.FindAll<string>(_s1.Split(
                        new char[] { ' ', '\t', ',' }), delegate (string s)
                        {
                            return !String.IsNullOrEmpty(s);
                        });
                string[] splitIntoWords2 = Array.FindAll<string>(_s2.Split(
                        new char[] { ' ', '\t', ',' }), delegate (string s)
                        {
                            return !String.IsNullOrEmpty(s);
                        });
                // [0] is loopID with length, [1] is clusterID, [2] is pdb and chainID, [6] is seq
                if (splitIntoWords1[0] != splitIntoWords2[0]) // loopIDs not equal
                {
                    if (splitIntoWords1[0].Substring(0, 1) != splitIntoWords2[0].Substring(0, 1))
                    { return splitIntoWords2[0].CompareTo(splitIntoWords1[0]); } // returns L before H
                    if (splitIntoWords1[0].Substring(1, 1) != splitIntoWords2[0].Substring(1, 1))
                    { return splitIntoWords1[0].CompareTo(splitIntoWords2[0]); } // returns L1 before L2
                                                                                 // if here, length is only difference
                    return splitIntoWords1[6].Length.CompareTo(splitIntoWords2[6].Length); //returns L1-10 before L1-11
                }
                if (splitIntoWords1[1] != splitIntoWords2[1]) // clusterIDs not equal
                { return splitIntoWords1[1].CompareTo(splitIntoWords2[1]); }
                return splitIntoWords1[2].CompareTo(splitIntoWords2[2]); // sorting pdbIDs in alphabetical order
            }
        }
        public class MyUniqueSequenceContainer
        {
            public MyUniqueSequenceContainer()
            {
                uniqueSeqObject = null;
                uniqueSeqCount = -1;
            }
            public MyUniqueSequenceContainer(string _seq)
            {
                uniqueSeqObject = (object)_seq;
                uniqueSeqCount = 1;
            }
            public MyUniqueSequenceContainer(string _seq, int _count)
            {
                uniqueSeqObject = (object)_seq;
                uniqueSeqCount = _count;
            }
            // memberVariables
            public object uniqueSeqObject = new object();
            public int uniqueSeqCount = new int();
        }
        public class MyUniqueSeqContainerSorter : IComparer<MyUniqueSequenceContainer>
        {
            public int Compare(MyUniqueSequenceContainer _s1, MyUniqueSequenceContainer _s2)
            {
                if (_s1.uniqueSeqCount != _s2.uniqueSeqCount)
                { return _s2.uniqueSeqCount.CompareTo(_s1.uniqueSeqCount); }
                else
                {
                    string seq1 = (string)_s1.uniqueSeqObject;
                    string seq2 = (string)_s2.uniqueSeqObject;
                    return seq1.CompareTo(seq2);
                }
            }
        }
        public class MyGeoSeqPair
        {
            public MyGeoSeqPair()
            {
                geoSeqObject = null;
                geoSeqCount = -1;
            }
            public MyGeoSeqPair(string _geoSeq, int _count)
            {
                geoSeqObject = (object)_geoSeq;
                geoSeqCount = _count;
            }
            // member variables
            public object geoSeqObject = new object();
            public int geoSeqCount = new int();
        }
        public class MyGeoSeqPairSorter : IComparer<MyGeoSeqPair>
        {
            public int Compare(MyGeoSeqPair _g1, MyGeoSeqPair _g2)
            {
                if (_g1.geoSeqCount != _g2.geoSeqCount)
                { return _g2.geoSeqCount.CompareTo(_g1.geoSeqCount); } // descending count order
                else
                {
                    string g1str = (string)_g1.geoSeqObject;
                    string g2str = (string)_g2.geoSeqObject;
                    return g1str.CompareTo(g2str); // if equal count, ascending alphabetical order
                }
            }
        }
        public class MyQualityObj
        {
            public MyQualityObj()
            {
                myPdbIDobj = null;
                myResolution = 999;
                myBfac = 999;
                myConfE = 999;
            }
            public MyQualityObj(string _pdbID, double _res, double _bfac, double _confE)
            {
                myPdbIDobj = (object)_pdbID;
                myResolution = _res;
                myBfac = _bfac;
                myConfE = _confE;
            }

            // memberVariables
            public object myPdbIDobj = new object();
            public double myResolution = new double();
            public double myBfac = new double();
            public double myConfE = new double();
        }

        public class MyQualitySorter : IComparer<MyQualityObj>
        {
            public int Compare(MyQualityObj _q1, MyQualityObj _q2)
            {
                // compare resolution, then bfactor, then confE
                double theTolerance = (double)0.05;
                if (Math.Abs(_q1.myResolution - _q2.myResolution) > theTolerance)
                { return _q1.myResolution.CompareTo(_q2.myResolution); } // sort in ascending order
                if (Math.Abs(_q1.myBfac - _q2.myBfac) > theTolerance)
                { return _q1.myBfac.CompareTo(_q2.myBfac); } // ditto
                if (Math.Abs(_q1.myConfE - _q2.myConfE) > theTolerance)
                { return _q1.myConfE.CompareTo(_q2.myConfE); } // ditto

                return 0;
            }
        }

        public class MyDihPDBpair
        {
            public MyDihPDBpair()
            {
                myPDBID = "XXXX";
                myDihedralValue = 999;
            }

            public MyDihPDBpair(string _pdbID, double _dihValue)
            {
                myPDBID = _pdbID;
                myDihedralValue = _dihValue;
            }
            // member variables
            public string myPDBID;
            public double myDihedralValue = new double();
        }

        public class MyDihSorter : IComparer<MyDihPDBpair>
        {
            public int Compare(MyDihPDBpair _p1, MyDihPDBpair _p2)
            {
                return _p1.myDihedralValue.CompareTo(_p2.myDihedralValue);
            }
        }
        public class MyClusterGroupSorter : IComparer<List<int>>
        {
            public int Compare(List<int> _a1, List<int> _a2)
            {
                int inta1 = _a1[0];
                int inta2 = _a2[0];
                return inta1.CompareTo(inta2);
            }
        }
        public class MyLoopKeySorter : IComparer<string>
        {
            public int Compare(string _l1, string _l2)
            {
                string loopType1 = _l1.Substring(0, 1);
                string loopType2 = _l2.Substring(0, 1);
                if (loopType1 != loopType2)
                { return _l2.CompareTo(_l1); } // should rank L before H
                string numeral1 = _l1.Substring(1, 1);
                string numeral2 = _l2.Substring(1, 1);
                if (numeral1 != numeral2)
                { return (Convert.ToInt32(numeral1)).CompareTo(Convert.ToInt32(numeral2)); }
                // next, find the length string and compare
                numeral1 = _l1.Substring(3, _l1.IndexOf("-", 3) - 3);
                numeral2 = _l2.Substring(3, _l2.IndexOf("-", 3) - 3);
                if (numeral1 != numeral2)
                { return (Convert.ToInt32(numeral1)).CompareTo(Convert.ToInt32(numeral2)); }
                return _l1.CompareTo(_l2); // this handles cis-trans hash            
            }
        }
        public class MyTurnDataObject
        {
            public MyTurnDataObject()
            {
                myTurnNameObj = (object)"empty";
                myTurnTurnIDObj = (object)"notDef";
                myClusterID = -1;
            }
            public MyTurnDataObject(string _name, string _loopKey, string _turnID, int _clID, ArrayList _dih)
            {
                myTurnNameObj = (object)_name;
                myLoopKeyObj = (object)_loopKey;
                myTurnTurnIDObj = (object)_turnID;
                myClusterID = _clID;
                myDihedrals = _dih;
            }
            // member variables
            public object myTurnNameObj = new object();
            public object myLoopKeyObj = new object();
            public object myTurnTurnIDObj = new object();
            public int myClusterID = new int();
            public ArrayList myDihedrals = new ArrayList();
        }
        #endregion


    }
}


