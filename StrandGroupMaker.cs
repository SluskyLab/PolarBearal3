using betaBarrelProgram.AtomParser;
using betaBarrelProgram.BarrelStructures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Numerics;

namespace betaBarrelProgram
{
    public class StrandGroupMaker //Rik
    {
        public string PdbName { get; set; }
        public List<GroupOfStrands> GroupOfGroup { get; set; }
        public List<Strand> BarrelStrands { get; set; }
        //public Chain myChain;

        public StrandGroupMaker(Protein newProt, List<Strand> AllStrands)
        {
            this.PdbName = newProt.PdbName;
            this.GroupOfGroup = new List<GroupOfStrands>();
            this.BarrelStrands = new List<Strand>();

            var UnassignedStrandSet = new List<Strand>(); //Create a set of unsorted strand with strand size greater than 3
            foreach (var strands in AllStrands)
            {
                if (strands.Residues.Count > 2)
                {
                    UnassignedStrandSet.Add(strands);
                }
            }

            var GroupedStrandSet = new GroupOfStrands();
            var FirstStrand = new OneStrand(UnassignedStrandSet[0]); // Selecting the first strand from the unsorted list of strands.
            UnassignedStrandSet.RemoveAt(0); //Remove the first strand from the unsorted list.
            GroupedStrandSet.StrandSet.Add(FirstStrand);

            do //Do this method till there are no more neighbors
            {

                var TotalNumberOfStrands = UnassignedStrandSet.Count;
                if ((UnassignedStrandSet.Count != 0) || (GroupedStrandSet.StrandSet.Count != 0))
                {
                    StrandSorter(UnassignedStrandSet, GroupedStrandSet.StrandSet); //Assign to different groups
                    if (TotalNumberOfStrands == (UnassignedStrandSet.Count/*-1*/))
                    {
                        Coordinate centroidOfGroup = Centroid(GroupedStrandSet); //Calculate the centroid of the group
                        GroupedStrandSet.Centroid = centroidOfGroup; //Add the centroid of the group
                        this.GroupOfGroup.Add(GroupedStrandSet); //Add the group of strands to the group of group.
                        FirstStrand = new OneStrand(UnassignedStrandSet[0]); // Selecting the first strand from the unsorted list of strands.
                        UnassignedStrandSet.RemoveAt(0); //Remove the first strand from the unsorted list.
                        GroupedStrandSet = new GroupOfStrands(); //Creat a new group
                        GroupedStrandSet.StrandSet.Add(FirstStrand); //Add the next first strand from the unassigned to the new group.                        
                    }
                }
                if ((UnassignedStrandSet.Count == 0) && (GroupedStrandSet.StrandSet.Count != 0)) //for the last strand set stuck in the GroupedStrandSet
                {
                    Coordinate centroidOfGroup = Centroid(GroupedStrandSet); //Calculate the centroid of the group
                    GroupedStrandSet.Centroid = centroidOfGroup; //Add the centroid of the group
                    this.GroupOfGroup.Add(GroupedStrandSet); //Add the group of strands to the group of group.
                    GroupedStrandSet = new GroupOfStrands(); //Creat a new group
                }
            } while ((UnassignedStrandSet.Count != 0) || (GroupedStrandSet.StrandSet.Count != 0));

            SharedFunctions.writePymolScriptForStrandGroups(GroupOfGroup, Global.OUTPUT_DIR, Global.DB_DIR, PdbName);

            #region Hardcoding the PDBs that fails to work on it's own
            var PDBList = new List<string>() {  "1O8V", "1IFC", "4UU3", "1GL4", "2HLV", "2YNK", "4GEY",
                                                "5GAQ", "6RB9", "3W9T", "5WC3",
                                                "6RHV", "2R73", "2RD7", "5AZO", "3ANZ",
                                                "4P1X", "5JZT",
                                                "1QWD", "3S26", "3HPE", "1R0U", "2JOZ", "4XMC", "2OOJ",
                                                "3KZA", "4R0B"};//#Hard Coding# Manually select strands
            bool ManualSelect = PDBList.Contains(PdbName);
            if (ManualSelect)
            {
                HardCode(GroupOfGroup);
            }
            #endregion

            foreach (var group in this.GroupOfGroup)
            {
                if ((group.StrandSet.Count >= 8) || (PdbName == "1VKY") || (PdbName == "1XQB"))
                {
                    if (!ManualSelect)
                    {
                        group.Cylinder = GetAxisCylinder(group);
                        NeighboringStrand(group);
                        CheckIfBarrel(group);
                    }
                    if (group.IsBarrel)
                    {
                        rearrangeBarrel(group);
                        AssignDirection(group);
                        BarrelAxis(group);
                        LogStrandData(group);
                        foreach (var strand in group.StrandSet)
                        {
                            BarrelStrands.Add(strand.StrandInTheGroup);
                        }
                        SharedFunctions.writePymolScriptForBarrelStrands(this.BarrelStrands, Global.OUTPUT_DIR, Global.DB_DIR, PdbName);
                        //SharedFunctions.WritePymolScriptForInOutStrands(this.BarrelStrands, Global.OUTPUT_DIR, Global.DB_DIR, PdbName, CCentroid, NCentroid);
                        break;
                    }
                }
            }


            //Console.WriteLine("");
        }

        private void HardCode(List<GroupOfStrands> groupOfGroup)
        {
            var pdbs = new Dictionary<string, int>() {  { "1O8V", 0 }, { "1IFC", 0 }, { "4UU3", 0 }, { "4GEY", 0 }, { "2HLV", 0 }, { "2R73", 0 }, 
                                                        { "2RD7", 0 }, { "1QWD", 0 }, { "3S26", 0 }, { "3HPE", 0 }, { "1R0U", 0 }, { "2JOZ", 0 }, 
                                                        { "4XMC", 0 }, { "2OOJ", 0 }, { "3KZA", 0 }, { "4R0B", 0 },
                                                        { "2YNK", 1 }, { "5GAQ", 1 }, { "1GL4", 2 }, { "6RB9", 3 }, { "5JZT", 3 }, { "3W9T", 3 }
            };//PDB name and the group number.

            if (new List<string>() { "2HLV", "2R73", "3S26", "3KZA", "4R0B" }.Contains(this.PdbName))
            {
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(8);//Removing problem strand
            }
            if (this.PdbName == "2RD7")
            {
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(8);//Removing problem strand
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(8);//Removing problem strand
            }
            if (this.PdbName == "1QWD")
            {
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(0);//Removing problem strand
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(0);//Removing problem strand
            }
            if (this.PdbName == "4XMC")
            {
                groupOfGroup[pdbs[this.PdbName]].StrandSet.RemoveAt(5);//Removing problem strand
            }


            if (new List<string>(pdbs.Keys).Contains(this.PdbName))
            {
                Console.WriteLine($"Manually changing 'IsBarrel' {this.PdbName}");
                var group = groupOfGroup[pdbs[this.PdbName]];
                group.IsBarrel = true;
            }
            else 
            {
                Console.WriteLine($"Manually didn't change anything for {this.PdbName}");
            }

        }

        public void AssignDirection(GroupOfStrands group)
        {
            for (int i = 0; i < group.StrandSet.Count; i++)
            {
                if (i == 0)
                {
                    group.StrandSet[0].DirectionAssigned = true;
                }
                else
                {
                    var strandOne = group.StrandSet[i - 1];
                    var strandTwo = group.StrandSet[i];
                    if (strandOne.DirectionAssigned = !strandTwo.DirectionAssigned)//check if either of the strands don't have a direction assigned
                    {
                        if (strandOne.DirectionAssigned == true) //Assign direction to the strand which doesn't have a direction yet
                        {
                            CheckIfParallel(strandTwo, strandOne);
                            strandTwo.DirectionAssigned = true;
                        }
                        else
                        {
                            CheckIfParallel(strandOne, strandTwo);
                            strandOne.DirectionAssigned = true;
                        }
                    }
                }

            }
        }

        private void BarrelAxis(GroupOfStrands group)
        {

            #region Selecting the ellipse based on which side of the plane the end residues are.
            List<Vector3> myCEllipse = new List<Vector3>();
            List<Vector3> myNEllipse = new List<Vector3>();
            Vector3 Ccentroid = new Vector3();
            Vector3 Ncentroid = new Vector3();
            int CCtr = 0;
            int NCtr = 0;
            int TCtr = 0;
            int ECtr = 0;

            foreach (var strand in group.StrandSet)
            {
                TCtr++;
                var strandStart = strand.StrandInTheGroup;
                var StartAtom = strandStart.Residues.Single(s => s.ResNum == strandStart.ResNumStart).BackboneCoords["CA"];
                var EndAtom = strandStart.Residues.Single(s => s.ResNum == strandStart.ResNumEnd).BackboneCoords["CA"];

                if (strand.DirectionAssigned)
                {
                    if (!strand.IsAntiParallel)
                    {
                        myCEllipse.Add(StartAtom);
                        Ccentroid += StartAtom;
                        CCtr++;
                        myNEllipse.Add(EndAtom);
                        Ncentroid += EndAtom;
                        NCtr++;
                    }
                    else
                    {
                        myCEllipse.Add(EndAtom);
                        Ccentroid += EndAtom;
                        CCtr++;
                        myNEllipse.Add(StartAtom);
                        Ncentroid += StartAtom;
                        NCtr++;
                    }
                }
                else
                {
                    ECtr++;
                    Console.WriteLine($"Direction not assigned for Strand {strand.StrandInTheGroup.betaStrandNum}!!");
                }
            }
            #endregion


            #region Defining the Ccentroid, Ncentroid and axis
            Ccentroid /= CCtr;
            Ncentroid /= NCtr;
            var Axis = Ncentroid - Ccentroid;
            //Console.WriteLine($"Ccentroid: {Ccentroid.X}, {Ccentroid.Y}, {Ccentroid.Z} and Ncentroid: {Ncentroid.X}, {Ncentroid.Y}, {Ncentroid.Z}");

            group.Ccentroid = Ccentroid;
            group.Ncentroid = Ncentroid;
            group.Axis = Axis;
            group.CellipseCoords = myCEllipse;
            group.NellipseCoords = myNEllipse;
            #endregion

            //foreach (var coordinate in myCEllipse)
            //{
            //    Console.WriteLine($"pseudoatom myCEllipse, pos=[{coordinate.X}, {coordinate.Y}, {coordinate.Z}]" +
            //        $"\nshow spheres, myCEllipse" +
            //        $"\ncolor red, myCEllipse " +
            //        $"");
            //}
            //foreach (var coordinate in myNEllipse)
            //{
            //    Console.WriteLine($"pseudoatom myNEllipse, pos=[{coordinate.X}, {coordinate.Y}, {coordinate.Z}]" +
            //        $"\nshow spheres, myNEllipse" +
            //        $"\ncolor blue, myNEllipse " +
            //        $"");
            //}


            #region Calculate InOut
            //List<Strand> strandList = new List<Strand>();
            //foreach (var strand in group.StrandSet)
            //{
            //    strandList.Add(strand.StrandInTheGroup);
            //}

            ////SharedFunctions.setInOut(strandList, Global.OUTPUT_DIR, this.PdbName, Axis, Ccentroid, Ncentroid);
            ////RYAN//SharedFunctions.setInOutMin(strandList, Global.OUTPUT_DIR, this.PdbName, Ccentroid, Ncentroid);
            //SharedFunctions.WritePymolScriptForInOutStrands(strandList, Global.OUTPUT_DIR, Global.DB_DIR, PdbName, Ccentroid, Ncentroid);
            ////SharedFunctions.setInOutMin(strandList, Global.OUTPUT_DIR, this.PdbName, Ccentroid, Ncentroid);
            #endregion
        }

        private void rearrangeBarrel(GroupOfStrands group)
        {
            Vector3 centroid = new Vector3((float)group.Centroid.X, (float)group.Centroid.Y, (float)group.Centroid.Z);
            List<Strand> strands = new List<Strand>();
            Vector3 firstCentroid = centroidCalculate(group.StrandSet[0].StrandInTheGroup);
            Vector3 reference = firstCentroid - centroid;

            foreach (var strand in group.StrandSet)
            {
                Vector3 strandCentroid = centroidCalculate(strand.StrandInTheGroup);

                Vector3 strandDirection = strandCentroid - centroid;
                strand.StrandInTheGroup.angle = SharedFunctions.AngleBetween(reference, strandDirection);

                strands.Add(strand.StrandInTheGroup);
            }

            Strand selectStrand = strands.Where(strand => (strand.angle > 60 && strand.angle < 130)).First();
            Vector3 selectStrandDirection = centroidCalculate(selectStrand) - centroid;
            Vector3 refCrossP = Vector3.Cross(reference, selectStrandDirection);

            var ctr = 0;
            foreach (var strand in group.StrandSet)
            {
                Vector3 strandCentroid = centroidCalculate(strand.StrandInTheGroup);
                Vector3 strandDirection = strandCentroid - centroid;
                Vector3 CrossP = Vector3.Cross(reference, strandDirection);
                var refAngle = SharedFunctions.AngleBetween(CrossP, refCrossP);
                if (refAngle > 90)
                {
                    strand.StrandInTheGroup.angle = 360 - strand.StrandInTheGroup.angle;
                }

                //Console.WriteLine($"pseudoatom strand{ctr}, pos=[{strandCentroid.X}, {strandCentroid.Y}, {strandCentroid.Z}]" +
                //    $"\nshow spheres, strand{ctr}" +
                //    $"\nlabel strand{ctr}, \"                    {strand.StrandInTheGroup.angle:000.##}\"" +
                //    $"");
                //ctr++;
            }

            //        Console.WriteLine($"pseudoatom centroid, pos=[{centroid.X}, {centroid.Y}, {centroid.Z}]" +
            //$"\nshow spheres, centroid");

            group.StrandSet.Sort((a, b) => a.StrandInTheGroup.angle.CompareTo(b.StrandInTheGroup.angle));

            //foreach (var item in strands)
            //{
            //    Console.WriteLine(item.angle);
            //}
            

            int i = 0;
            foreach (var strand in group.StrandSet)
            {
                //Console.WriteLine($"{item.betaStrandNum} angle is {item.angle}");
                strand.StrandInTheGroup.StrandNum = i;
                i++;
            }

            Vector3 centroidCalculate(Strand Strand)
            {
                var caList = Strand.Select(residue => residue.BackboneCoords["CA"]).ToList();
                Vector3 strandCentroid = new Vector3
                {
                    X = caList.Sum(ca => ca.X) / caList.Count,
                    Y = caList.Sum(ca => ca.Y) / caList.Count,
                    Z = caList.Sum(ca => ca.Z) / caList.Count
                };
                return strandCentroid;
            }
            //return strands;
        }

        private void CheckIfBarrel(GroupOfStrands group)
        {


            int NeighborCount = 0;
            foreach (var oneStrand in group.StrandSet)
            {
                NeighborCount += oneStrand.LeftNeighbors.Count;
                NeighborCount += oneStrand.RightNeighbors.Count;
            }
            if (NeighborCount == (2 * group.StrandSet.Count)) //Closed beta barrels should have 2n neighbors compared to 2n-2 for open beta sheet
            {
                group.IsBarrel = true;
            }
            if ((NeighborCount == ((2 * group.StrandSet.Count) - 2)) && (group.StrandSet.Count > 4))
            {
                {
                    foreach (var oneStrand in group.StrandSet)
                    {
                        if (oneStrand.LeftNeighbors.Count == 0 || oneStrand.RightNeighbors.Count == 0)
                        {
                            OneStrand strand1 = oneStrand;
                            foreach (var oneStrand2 in group.StrandSet)
                            {
                                if ((oneStrand2 != strand1) && (oneStrand2.LeftNeighbors.Count == 0 || oneStrand2.RightNeighbors.Count == 0))
                                {
                                    OneStrand strand2 = oneStrand2;
                                    group.IsBarrel = CheckIfNear(strand1, strand2);
                                    break;
                                }
                            }
                            break;
                        }
                    }

                }
            }
        }

        private bool CheckIfNear(OneStrand strand1, OneStrand strand2)
        {
            bool flag = false;
            foreach (var Res1 in strand1.StrandInTheGroup.Residues)
            {
                foreach (var Res2 in strand2.StrandInTheGroup.Residues)
                {
                    var distance = Res1.BackboneCoords["CA"] - Res2.BackboneCoords["CA"];
                    if (distance.Length() < 11)
                    {
                        flag = true;
                    }
                }
            }
            return flag;
        }

        private void StrandSorter(List<Strand> UnassignedStrandSet, List<OneStrand> GroupedStrandSet)
        {
            for (int i = 0; i < UnassignedStrandSet.Count; i++)
            {

                for (int j = 0; j < GroupedStrandSet.Count; j++)
                {
                    var SortedStrand = GroupedStrandSet[j];
                    int Hbond = 0;
                    if (Hbond < 1)
                    {
                        for (int r = 0; r < UnassignedStrandSet[i].Residues.Count; r++)
                        {
                            Res Residue1 = UnassignedStrandSet[i].Residues[r];
                            //Console.WriteLine($"Residue {r} in unassigned strand {i}, total residues in strand: {UnassignedStrandSet[i].Residues.Count}");

                            for (int p = 0; p < GroupedStrandSet[j].StrandInTheGroup.Residues.Count; p++)
                            {
                                Res Residue2 = GroupedStrandSet[j].StrandInTheGroup.Residues[p];
                                //Console.WriteLine($"Residue {p} in sorted strand {j}, total resideue in strand: {GroupedStrandSet[j].StrandInTheGroup.Residues.Count}");


                                if (CheckHBond(Residue1, Residue2, UnassignedStrandSet[i], GroupedStrandSet[j].StrandInTheGroup)) //Pass the two residue to the function CheckHBond to see if they form Hbond.
                                {
                                    Hbond++;
                                    if (Hbond >= 1)
                                    {
                                        break;
                                    }
                                }
                            }
                            if (Hbond >= 1)
                            {
                                break;
                            }
                        }
                    }
                    if (Hbond >= 1)
                    {

                        //Console.WriteLine($"\n\nUnassigned count = {UnassignedStrandSet.Count}");
                        //SortedStrand.AllNeighboringStrands.Add(UnassignedStrandSet[i]); //Add the neighboring strand as the neighbor of the sorted strand.
                        var TempStrand = new OneStrand(UnassignedStrandSet[i]); //Create a new groupedstrand with the neighbor strand.
                        GroupedStrandSet.Add(TempStrand); //Add the neighboring grouped strand to the group.
                        UnassignedStrandSet.RemoveAt(i); //Remove that strand from the unassigned.
                        i = 0;
                        j = 0; //Reset counters as the the count has modified.
                        if ((UnassignedStrandSet.Count == 0))
                        {
                            break;
                        }
                        //if (UnassignedStrandSet.Count <= i)
                        //{
                        //    i = -1;
                        //}
                        //if (GroupedStrandSet.Count <= j)
                        //{
                        //    j = -1;
                        //}

                        #region CodeTest
                        //Console.WriteLine($"i = {i}, j = {j}");
                        //Console.WriteLine($" tempstrand = {TempStrand.StrandInTheGroup.ResNumStart}");
                        //Console.WriteLine($" unsortedStrand = {UnsortedStrand.ResNumStart}");
                        //Console.WriteLine($" unassigned at i = {UnassignedStrandSet[i].ResNumStart}");

                        //Console.WriteLine("\nUnassignedStrandSet\n");
                        //for (int a = 0; a < UnassignedStrandSet.Count; a++)
                        //{
                        //    Console.WriteLine($" Pos{a} : {UnassignedStrandSet[a].ResNumStart}");
                        //}
                        //Console.WriteLine("\nGroupedStrandSet\n");
                        //for (int a = 0; a < GroupedStrandSet.Count; a++)
                        //{
                        //    Console.WriteLine($" Pos{a} : {GroupedStrandSet[a].StrandInTheGroup.ResNumStart}");
                        //}
                        #endregion
                    }
                }
                if (UnassignedStrandSet.Count == 0)
                {
                    break;
                }
            }
        }

        public void NeighboringStrand(GroupOfStrands groupOfStrands)
        {
            Coordinate centroidOfGroup = groupOfStrands.Centroid;
            int k = 0;

            groupOfStrands.StrandSet[0].DirectionAssigned = true; //Assigning the direction of first strand as parrallel
            for (int i = 0; i < groupOfStrands.StrandSet.Count; i++) //Iterate over two strand and check if they are neigbors based on Hbonds
            {
                OneStrand strandOne = groupOfStrands.StrandSet[i];
                strandOne.LeftNeighbors = new List<Strand>();
                strandOne.RightNeighbors = new List<Strand>();

                for (int j = 0; j < groupOfStrands.StrandSet.Count; j++)
                {
                    OneStrand strandTwo = groupOfStrands.StrandSet[j];

                    k++;
                    int Hbond = 0;
                    var resList1 = new List<Res>();
                    var resList2 = new List<Res>();
                    //Console.WriteLine($"StrandOne: {strandOne.StrandInTheGroup.betaStrandNum}, strandTwo: {strandTwo.StrandInTheGroup.betaStrandNum}");
                    if (strandOne != strandTwo)
                    {
                        for (int p = 0; p < strandOne.StrandInTheGroup.Residues.Count; p++)
                        {
                            Res ResOne = strandOne.StrandInTheGroup.Residues[p];
                            for (int r = 0; r < strandTwo.StrandInTheGroup.Residues.Count; r++)
                            {
                                Res ResTwo = strandTwo.StrandInTheGroup.Residues[r];
                                //Console.WriteLine($"i:{i}, j:{j},p:{p},r:{r}");
                                if (CheckHBond(ResOne, ResTwo, strandOne.StrandInTheGroup, strandTwo.StrandInTheGroup))
                                {
                                    Hbond++;
                                    if (resList1.Count == 0)
                                    {
                                        resList1.Add(ResOne);
                                    }
                                    else if (resList1.Any(item => item.SeqID != ResOne.SeqID)) //make sure the residue is already not in the list  
                                    {
                                        if (((resList1.Min(item => (item.BackboneCoords["CA"] - ResOne.BackboneCoords["CA"]).Length())) < 10)) //The residue is not too far away from residues in the list
                                        {
                                            resList1.Add(ResOne);
                                        }
                                        else
                                        {
                                            var dist = (resList1.Min(item => (item.BackboneCoords["CA"] - ResOne.BackboneCoords["CA"]).Length()));
                                            //Console.WriteLine($"ResOne Distance: {dist}");
                                        }
                                    }
                                    if (resList2.Count == 0)
                                    {
                                        resList2.Add(ResTwo);
                                    }
                                    else if (resList2.Any(item => item.SeqID != ResTwo.SeqID))//make sure the residue is already not in the list
                                    {
                                        if (((resList2.Min(item => (item.BackboneCoords["CA"] - ResTwo.BackboneCoords["CA"]).Length())) < 10)) //The residue is not too far away from residues in the list
                                        {
                                            resList2.Add(ResTwo);
                                        }
                                        else
                                        {
                                            var dist = (resList2.Min(item => (item.BackboneCoords["CA"] - ResTwo.BackboneCoords["CA"]).Length()));
                                            //Console.WriteLine($"ResTwo Distance: {dist}");
                                        }
                                    }
                                }
                                //if (Hbond > 0)
                                //{
                                //    break;
                                //}
                            }
                            //if (Hbond > 0)
                            //{
                            //    break;
                            //}
                        }
                    }
                    if (Hbond != 0) //if the the strands are neighbors then check if the second strand is left or right of the first one.
                    {
                        //using (StreamWriter file = new StreamWriter(fileLocation, append: true))
                        //{
                        //    string newLine = this.PdbName + "\t" + strandOne.StrandInTheGroup.betaStrandNum + "\t" + strandTwo.StrandInTheGroup.betaStrandNum + "\t" + Hbond; // "FirstAtom" + "\t" + "LastAtom" + "\t" + "StrandTwoAtom";
                        //    file.Write(newLine);
                        //}
                        if (strandOne.StrandInTheGroup.betaStrandNum == 7)
                        {
                            //Console.ReadKey();
                        }
                        //Double detTest = CheckLeftOrRightStrandResList(strandOne.StrandInTheGroup, strandTwo.StrandInTheGroup, resList1, resList2, centroidOfGroup);
                        Double determinant = 0;
                        if (!groupOfStrands.Cylinder.Collision)
                        {
                            determinant = CheckLeftOrRightCylinder(resList1, resList2, groupOfStrands); //Use the method to determine left or right 
                        }
                        else
                        {
                            determinant = CheckLeftOrRightResList(resList1, resList2, centroidOfGroup); //Use the method to determine left or right
                        }

                        if (determinant < 0) //If determinant is -ve then the strand is on the left.
                        {
                            strandOne.LeftNeighbors.Add(strandTwo.StrandInTheGroup);
                        }
                        else
                        {
                            strandOne.RightNeighbors.Add(strandTwo.StrandInTheGroup);
                        }
                        //if (strandOne.DirectionAssigned = !strandTwo.DirectionAssigned)//check if either of the strands don't have a direction assigned
                        //{
                        //    if (strandOne.DirectionAssigned == true) //Assign direction to the strand which doesn't have a direction yet
                        //    {
                        //        CheckIfParallel(strandTwo, strandOne);
                        //        strandTwo.DirectionAssigned = true;
                        //    }
                        //    else
                        //    {
                        //        CheckIfParallel(strandOne, strandTwo);
                        //        strandOne.DirectionAssigned = true;
                        //    }
                        //} //12/18/2020 - Removing strand assignment here. will implement it elsewhere.

                    }

                }
                //Console.WriteLine($"Strand One Left Count: {strandOne.LeftNeighbors.Count}, Strand One Right Count: {strandOne.RightNeighbors.Count} ");
                if ((strandOne.LeftNeighbors.Count > 1) || (strandOne.RightNeighbors.Count > 1))
                {
                    int ResCutOff = 2; //Ignore neighboring strands with 2 or less residues.
                    //check if any of the neighboring strand is less than 4 Residues
                    if (strandOne.LeftNeighbors.Count > 1)
                    {
                        foreach (var item in strandOne.LeftNeighbors)
                        {
                            if (item.NumOfRes <= ResCutOff)
                            {
                                groupOfStrands.StrandSet.Remove(groupOfStrands.StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == item.betaStrandNum));
                                strandOne.LeftNeighbors.Remove(item);
                                i = -1;
                                break;
                            }
                        }
                    }
                    if (strandOne.RightNeighbors.Count > 1)
                    {
                        foreach (var item in strandOne.RightNeighbors)
                        {
                            if (item.NumOfRes <= ResCutOff)
                            {
                                groupOfStrands.StrandSet.Remove(groupOfStrands.StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == item.betaStrandNum));
                                strandOne.RightNeighbors.Remove(item);
                                i = -1;
                                break;
                            }
                        }
                    }
                    if ((strandOne.LeftNeighbors.Count > 1) || (strandOne.RightNeighbors.Count > 1))
                    {
                        #region Uncomment this region and add breakpoint to see check strand selection
                        var ListOfStrand = new List<Strand>();
                        foreach (var oneStrand in groupOfStrands.StrandSet)
                        {
                            ListOfStrand.Add(oneStrand.StrandInTheGroup);
                        }
                        //SharedFunctions.writePymolScriptForBarrelStrands(ListOfStrand, Global.OUTPUT_DIR, Global.DB_DIR, PdbName);
                        #endregion
                        Concatenate(groupOfStrands);
                        i = -1;
                    }
                }
            }
        }

        private void CheckIfParallel(OneStrand strand1, OneStrand strand2)//first strand needs a direction assigned
        {
            var resiOne = strand1.StrandInTheGroup.Residues;
            var resiTwo = strand2.StrandInTheGroup.Residues;

            var resFirstStrand1 = resiOne[0].BackboneCoords["CA"];
            var resFirstStrand2 = resiTwo[0].BackboneCoords["CA"];
            var resLastStrand1 = resiOne[resiOne.Count - 1].BackboneCoords["CA"];
            var resLastStrand2 = resiTwo[resiTwo.Count - 1].BackboneCoords["CA"];

            var vectorStrand1 = resFirstStrand1 - resLastStrand1;
            var vectorStrand2 = resFirstStrand2 - resLastStrand2;

            var angleBetween = SharedFunctions.AngleBetween(vectorStrand1, vectorStrand2);

            //Console.WriteLine(angleBetween);

            if ((strand2.IsAntiParallel == false && angleBetween >= 90) || (strand2.IsAntiParallel == true && angleBetween < 90))
            {
                strand1.IsAntiParallel = true;
            }
        }

        // Calculate the centroid of a group of Strands based on the CA of each residue.  
        private Coordinate Centroid(GroupOfStrands groupOfStrands)
        {
            // !!!ATENTION!!! MAybe change this later to read 
            int i = 0;
            double x = 0;
            double y = 0;
            double z = 0;
            foreach (OneStrand strand in groupOfStrands.StrandSet)
            {
                foreach (Res residue in strand.StrandInTheGroup.Residues)
                {
                    foreach (Atom atom in residue)
                    {
                        if (atom.AtomName == "CA")
                        {
                            i++;
                            x = x + atom.Coords.X;
                            y = y + atom.Coords.Y;
                            z = z + atom.Coords.Z;
                        }
                    }
                }
            }
            Coordinate centroid = new Coordinate((x / i), (y / i), (z / i));
            return centroid;
        }

        private void Concatenate(GroupOfStrands groupOfStrands)
        {
            var StrandSet = groupOfStrands.StrandSet;
            for (int i = 0; i < StrandSet.Count; i++)
            {

                if (StrandSet[i].LeftNeighbors.Count > 1) //For any strand that has more than 1 neighbor.
                {
                    //moveNeighbor(StrandSet[i].LeftNeighbors);
                    joinNeighbors(StrandSet[i].LeftNeighbors);
                    if (StrandSet[i].LeftNeighbors.Count < 2 && StrandSet[i].RightNeighbors.Count < 2)
                    {
                        break;
                    }
                }
                if (StrandSet[i].RightNeighbors.Count > 1) //For any strand that has more than 1 neighbor.
                {
                    //moveNeighbor(StrandSet[i].RightNeighbors);
                    joinNeighbors(StrandSet[i].RightNeighbors);
                    if (StrandSet[i].LeftNeighbors.Count < 2 && StrandSet[i].RightNeighbors.Count < 2)
                    {
                        break;
                    }
                }



                //8/11/2020 - if there are more than 2 neighbors, this fails and causes an infinite loop
                //void moveNeighbor(List<Strand> Neighbors)
                //{
                //    for (int j = (Neighbors.Count - 1); j == 1; j--)
                //    {
                //        var StrandNum1 = Neighbors[j].betaStrandNum;
                //        var StrandNum2 = Neighbors[j - 1].betaStrandNum;
                //        var Strand1 = StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == StrandNum1).StrandInTheGroup;
                //        var Strand2 = StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == StrandNum2).StrandInTheGroup;
                //        Strand1.Residues.AddRange(Strand2.Residues);

                //        StartEnd(Strand1, Strand2);
                //        //Strand1.ResNumEnd = Strand2.ResNumEnd;

                //        Strand1.NumOfRes = Strand1.Residues.Count;
                //        StrandSet.Remove(StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == StrandNum2));

                //    }
                //}


                //8/11/2020 - Adding this new method to deal with more than 2 neighbors
                void joinNeighbors(List<Strand> Neighbors)
                {
                    do
                    {
                        Strand firstStrand = findfarthestNeighbor(Neighbors);
                        Strand nearestStrand = findNearestNeighbor(firstStrand, Neighbors);
                        firstStrand.Residues.AddRange(nearestStrand.Residues);
                        StartEnd(firstStrand, nearestStrand);
                        //firstStrand.Residues.Sort((x, y) => x.SeqID.CompareTo(y.SeqID)); //8/23/20 - Removed, as it should be sort by physical  position. //This sorts residue by sequence ID 
                        var firstRes = firstStrand.Residues.Single(s => s.ResNum == firstStrand.ResNumStart);
                        firstStrand.Residues.Sort((x, y) =>
                        {
                            var d1 = (firstRes.BackboneCoords["CA"] - x.BackboneCoords["CA"]).Length();
                            var d2 = (firstRes.BackboneCoords["CA"] - y.BackboneCoords["CA"]).Length();
                            //Console.WriteLine($"d1 = {d1}, d2= {d2}");
                            return d1.CompareTo(d2);
                        });//sort residues based on there distance from the first strand
                        StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == firstStrand.betaStrandNum).Concatenated = true;
                        StrandSet.Remove(StrandSet.Single(s => s.StrandInTheGroup.betaStrandNum == nearestStrand.betaStrandNum));
                        Neighbors.Remove(Neighbors.Single(s => s.betaStrandNum == nearestStrand.betaStrandNum));
                    } while (Neighbors.Count > 1);
                }

                //method to find the strand farthest from it's neighbor
                Strand findfarthestNeighbor(List<Strand> neighbors)
                {
                    double maxOfMaxLength = 0;
                    Strand maxStrand = neighbors[0];
                    foreach (var strand1 in neighbors)
                    {
                        foreach (var strand2 in neighbors)
                        {
                            if (strand1 != strand2)
                            {
                                var listOfDist = CompareStrands(strand1, strand2);
                                var maxLength = listOfDist.Max(r => r.Length());
                                if (maxLength > maxOfMaxLength)//change the strand and longest length if that is found
                                {
                                    maxOfMaxLength = maxLength;
                                    maxStrand = strand1;
                                }
                            }
                        }
                    }
                    return maxStrand;
                }

                //method to find the nearest strand to the first strand
                Strand findNearestNeighbor(Strand firstStrand, List<Strand> neighbors)
                {
                    double minOfMinLength = 0;
                    Strand minStrand = neighbors[neighbors.Count - 1];
                    foreach (var strand in neighbors)
                    {
                        if (strand != firstStrand)
                        {
                            var listOfDist = CompareStrands(firstStrand, strand);
                            var minLength = listOfDist.Min(r => r.Length());
                            if (minLength > minOfMinLength)//change the strand and longest length if that is found
                            {
                                minOfMinLength = minLength;
                                minStrand = strand;
                            }
                        }
                    }
                    return minStrand;
                }

                //method to get the list of distance between the first and last residue of two strands
                List<Vector3> CompareStrands(Strand strand1, Strand strand2)
                {
                    var PointOneA = strand1.Residues.Single(s => s.ResNum == strand1.ResNumStart).BackboneCoords["CA"];
                    var PointOneB = strand1.Residues.Single(s => s.ResNum == strand1.ResNumEnd).BackboneCoords["CA"];
                    var PointTwoA = strand2.Residues.Single(s => s.ResNum == strand2.ResNumStart).BackboneCoords["CA"];
                    var PointTwoB = strand2.Residues.Single(s => s.ResNum == strand2.ResNumEnd).BackboneCoords["CA"];

                    var AA = PointOneA - PointTwoA;
                    var AB = PointOneA - PointTwoB;
                    var BA = PointOneB - PointTwoA;
                    var BB = PointOneB - PointTwoB;

                    var listOfDist = new List<Vector3> { AA, AB, BA, BB };
                    return listOfDist;
                }


                void StartEnd(Strand strand1, Strand strand2)
                {
                    var listOfDist = CompareStrands(strand1, strand2);
                    double maxLength = listOfDist.Max(r => r.Length());

                    //Console.WriteLine($"Start: {strand1.Residues.Single(s => s.ResNum == strand1.ResNumStart).SeqID} & End: {strand1.Residues.Single(s => s.ResNum == strand1.ResNumEnd).SeqID}");

                    if (listOfDist[0].Length() == maxLength)
                    {
                        strand1.ResNumEnd = strand2.ResNumStart;
                    }
                    else if (listOfDist[1].Length() == maxLength)
                    {
                        strand1.ResNumEnd = strand2.ResNumEnd;
                    }
                    else if (listOfDist[2].Length() == maxLength)
                    {
                        strand1.ResNumStart = strand2.ResNumStart;
                    }
                    else if (listOfDist[3].Length() == maxLength)
                    {
                        strand1.ResNumStart = strand2.ResNumEnd;
                    }

                    //Console.WriteLine($"Start: {strand1.Residues.Single(s => s.ResNum == strand1.ResNumStart).SeqID} & End: {strand1.Residues.Single(s => s.ResNum == strand1.ResNumEnd).SeqID}");
                }
            }

        }

        private double CheckLeftOrRightResList(List<Res> ResList1, List<Res> ResList2, Coordinate centroid)
        {

            Vector3 FirstAtom = new Vector3(); //A
            Vector3 LastAtom = new Vector3(); //B
            Vector3 Centroid = new Vector3((float)centroid.X, (float)centroid.Y, (float)centroid.Z); //C
            Vector3 StrandTwoAtom = new Vector3(); //X

            //Console.WriteLine($"ResList1 Count : {ResList1.Count}, ResList Count2 : {ResList2.Count}");

            ////Display the two list
            //Console.WriteLine("ResList1:");
            //foreach (var res in ResList1)
            //{
            //    Console.Write($"{res.SeqID} ");
            //}
            //Console.WriteLine("");
            //Console.WriteLine("ResList2:");
            //foreach (var res in ResList2)
            //{
            //    Console.Write($"{res.SeqID} ");
            //}
            //Console.WriteLine("");

            var medianResList1 = (ResList1.Count) / 2;
            var medianResList2 = (ResList2.Count) / 2;

            if (ResList1.Count >= 3)
            {
                FirstAtom = ResList1[medianResList1 - 1].BackboneCoords["CA"];
                LastAtom = ResList1[medianResList1 + 1].BackboneCoords["CA"];
                StrandTwoAtom = ResList2[medianResList2].BackboneCoords["CA"];

                //string fileLocation = Global.OUTPUT_DIR + "NeighboringStrandLog.txt";
                //using (StreamWriter file = new StreamWriter(fileLocation, append: true))
                //{
                //    string newLine = "\t" + ResList1[medianResList1 - 1].SeqID + "\t" + ResList1[medianResList1 + 1].SeqID + "\t" + ResList2[medianResList2].SeqID;
                //    file.Write(newLine);
                //}
            }
            else
            {
                if (ResList1.Count == 2)
                {
                    FirstAtom = ResList1[0].BackboneCoords["CA"];
                    LastAtom = ResList1[1].BackboneCoords["CA"];
                    StrandTwoAtom = ResList2[0].BackboneCoords["CA"];

                    //string fileLocation = Global.OUTPUT_DIR + "NeighboringStrandLog.txt";
                    //using (StreamWriter file = new StreamWriter(fileLocation, append: true))
                    //{
                    //    string newLine = "\t" + ResList1[0].SeqID + "\t" + ResList1[1].SeqID + "\t" + ResList2[0].SeqID;
                    //    file.Write(newLine);
                    //}
                }
                else
                {
                    FirstAtom = ResList1[0].BackboneCoords["N"];
                    LastAtom = ResList1[0].BackboneCoords["C"];
                    StrandTwoAtom = ResList2[0].BackboneCoords["CA"];

                    //string fileLocation = Global.OUTPUT_DIR + "NeighboringStrandLog.txt";
                    //using (StreamWriter file = new StreamWriter(fileLocation, append: true))
                    //{
                    //    string newLine = "\t" + ResList1[0].SeqID + "\t" + ResList1[0].SeqID + "\t" + ResList2[0].SeqID;
                    //    file.Write(newLine);
                    //}
                }
            }

            var determinantPoint = Determinant(FirstAtom, LastAtom, StrandTwoAtom, Centroid);
            //Console.WriteLine(determinantPoint);
            return determinantPoint;

        }

        //Method to return a positive or negative value based on which side of the plane the MEDIAN strand residue is at.
        private double CheckLeftOrRightCylinder(List<Res> ResList1, List<Res> ResList2, GroupOfStrands group)
        {
            Vector3 FirstAtom, LastAtom, StrandTwoAtom, Reference;

            var medianResList1 = (ResList1.Count) / 2;
            var medianResList2 = (ResList2.Count) / 2;

            if (ResList1.Count >= 3) //When there are more than three atoms, use spaced out CA atoms.
            {
                FirstAtom = ResList1[medianResList1 - 1].BackboneCoords["CA"];
                LastAtom = ResList1[medianResList1 + 1].BackboneCoords["CA"];
                StrandTwoAtom = ResList2[medianResList2].BackboneCoords["CA"];
            }
            else
            {
                if (ResList1.Count == 2) //When there are only two atoms, use those two CA atoms.
                {
                    FirstAtom = ResList1[0].BackboneCoords["CA"];
                    LastAtom = ResList1[1].BackboneCoords["CA"];
                    StrandTwoAtom = ResList2[0].BackboneCoords["CA"];
                }
                else //When there is only one atom, use the N and C atoms.
                {
                    FirstAtom = ResList1[0].BackboneCoords["N"];
                    LastAtom = ResList1[0].BackboneCoords["C"];
                    StrandTwoAtom = ResList2[0].BackboneCoords["CA"];
                }
            }

            //Calculate nearest reference point on the cylinder axis
            Cylinder cylinder = group.Cylinder;
            Vector3 p1 = cylinder.pointOne;
            Vector3 p2 = cylinder.pointTwo;
            Vector3 q = FirstAtom;
            Vector3 u = p2 - p1;
            Vector3 pq = q - p1;
            Vector3 w2 = pq - Vector3.Multiply(u, Vector3.Dot(pq, u) / u.LengthSquared());
            Reference = q - w2;

            //Output a file with list of Pseudoatom with reference
            //SharedFunctions.writePymolScriptReference(group, FirstAtom, LastAtom, StrandTwoAtom, Reference, Global.OUTPUT_DIR, Global.DB_DIR, PdbName);

            var determinantPoint = Determinant(FirstAtom, LastAtom, StrandTwoAtom, Reference);
            //Console.WriteLine(determinantPoint);
            return determinantPoint;

        }

        private bool CheckHBond(Res residue1, Res residue2, Strand strand1, Strand strand2)
        {
            double d;
            double d1 = 0;
            double d2 = 0;
            double MinDist = 3.5;
            //foreach (Res residue1 in AllStrands[10].Residues)
            //{
            //    foreach (Res residue2 in AllStrands[11].Residues)
            //    {
            //int n = 0;
            //int m = 0;
            foreach (Atom atom1 in residue1.Atoms)
            {
                //n++;

                foreach (Atom atom2 in residue2.Atoms)
                {
                    //m++;
                    //Console.WriteLine($"Loop 1 : {n}, Loop : {m}");
                    var A1 = atom1.AtomName;
                    var A2 = atom2.AtomName;
                    //Console.WriteLine($"{A1} of {residue1.SeqID}{residue1.ThreeLetCode} & {A2} of {residue2.SeqID}{residue2.ThreeLetCode}");
                    if (((A1 == "O") && (A2 == "N")))
                    {
                        d1 = (atom1.Coords - atom2.Hydrogen).Length();

                        //if (d1 < MinDist)
                        //{
                        //    //break;
                        //    Console.WriteLine($"{A1} of {residue1.SeqID}{residue1.ThreeLetCode} of strand {strand1.betaStrandNum} forms bond with {A2} of {residue2.SeqID}{residue2.ThreeLetCode} of strand {strand2.betaStrandNum} at {d1}");
                        //}
                        ////else
                        ////{
                        ////    Console.WriteLine($"ELSE: {A1} of {residue1.SeqID}{residue1.ThreeLetCode} of strand {strand1.betaStrandNum} forms bond with {A2} of {residue2.SeqID}{residue2.ThreeLetCode} of strand {strand2.betaStrandNum} at {d1}");
                        ////}

                    }
                    if ((A1 == "N") && (A2 == "O"))
                    {
                        d2 = (atom2.Coords - atom1.Hydrogen).Length();


                        //if (d2 < MinDist)
                        //{
                        //    //break;
                        //    Console.WriteLine($"{A1} of {residue1.SeqID}{residue1.ThreeLetCode} of strand {strand1.betaStrandNum} forms bond with {A2} of {residue2.SeqID}{residue2.ThreeLetCode} of strand {strand2.betaStrandNum} at {d2}");
                        //}
                        ////else
                        ////{
                        ////    Console.WriteLine($"ELSE: {A1} of {residue1.SeqID}{residue1.ThreeLetCode} of strand {strand1.betaStrandNum} forms bond with {A2} of {residue2.SeqID}{residue2.ThreeLetCode} of strand {strand2.betaStrandNum} at {d2}");
                        ////}

                    }


                    //if (false)//((A1 == "O")&& (A2 == "N"))
                    //{
                    //    double d = (atom1.Coords - atom2.Hydrogen).Length;
                    //    double w = SharedFunctions.calculateHydrogenBond(atom2, atom1, d);
                    //    Console.WriteLine($"{A1} of {residue1.ThreeLetCode} of strand {AllStrands[0].ResNumStart} forms bond with {A2} of {residue2.ThreeLetCode} of strand {AllStrands[1].ResNumStart} \nw = {w}\n");
                    //}
                    //if (false)//((A1 == "N") && (A2 == "O"))
                    //{
                    //    double d = (atom2.Coords - atom1.Hydrogen).Length;
                    //    double w = SharedFunctions.calculateHydrogenBond(atom1, atom2, d);
                    //    Console.WriteLine($"{A1} of {residue1.ThreeLetCode} of strand {AllStrands[0].ResNumStart} forms bond with {A2} of {residue2.ThreeLetCode} of strand {AllStrands[1].ResNumStart} \nw = {w}\n");
                    //}


                }

            }

            //    }

            //}
            if (d1 <= d2)
            {
                d = d1;
            }
            else
            {
                d = d2;
            }
            return (d < MinDist);

        }

        // Determinant calculation for plane containing A, B, C and checking position of X.
        private double Determinant(Vector3 A, Vector3 B, Vector3 X, Vector3 C)
        {
            Vector3 Asub = Vector3.Subtract(A, C); //(A' = A - C)
            Vector3 Bsub = Vector3.Subtract(B, C); //(B' = B - C)
            Vector3 Xsub = Vector3.Subtract(X, C); //(X' = X - C)

            var a = Asub.X;
            var b = Asub.Y;
            var c = Asub.Z;
            var d = Bsub.X;
            var e = Bsub.Y;
            var f = Bsub.Z;
            var g = Xsub.X;
            var h = Xsub.Y;
            var i = Xsub.Z;

            var determinantPoint = a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;

            return determinantPoint;
        }

        public Cylinder GetAxisCylinder(GroupOfStrands group)
        {
            //Make a list of CA atom to work with.
            var caList = group.StrandSet.SelectMany(oneStrand => oneStrand.StrandInTheGroup.Select(residue => residue.BackboneCoords["CA"])).ToList();

            //Calculate centroid of the CA atom list.
            Vector3 caCentroid = new Vector3
            {
                X = caList.Sum(ca => ca.X) / caList.Count,
                Y = caList.Sum(ca => ca.Y) / caList.Count,
                Z = caList.Sum(ca => ca.Z) / caList.Count
            };

            //Move all the CA Atom coordinates around the origin from the original centroid
            for (int i = 0; i < caList.Count; i++)
            {
                caList[i] = caList[i] - caCentroid;
            }

            //Create a cylinder
            Cylinder cylinder = new Cylinder
            {
                pointOne = new Vector3((caList.Max(ca => ca.Length()) + 2), 0, 0),
                pointTwo = new Vector3(-(caList.Max(ca => ca.Length()) + 2), 0, 0),
                radius = 2.5
            };
            cylinder.moment = Vector3.Cross(cylinder.pointOne, cylinder.pointTwo);

            //outputcylinder();
            //Console.WriteLine("");

            #region Fourth attempt

            //    bool firstPass = false;
            //secondPass:
            //Create Variable
            double degree = 1;
            double radiusChange = 0.5;
            double totalAngle = 0;
            Dictionary<double, bool> xAxis, yAxis, zAxis;
            double xAngle, yAngle, zAngle;

            //outputcylinder();
            //Console.WriteLine("");

            //Second rotation check in Y axis
            yAxis = CollisionAngles(new Vector3(0, (float)degree, 0)); //creates a dictionary of angles and if the angles causes collision 
            totalAngle = yAxis.Keys.Max();
            Rotate(cylinder, new Vector3(0, -(float)totalAngle, 0)); //rotate back to original position
            List<double> listY = yAxis.Where(kvp => kvp.Value == false).Select(kvp => kvp.Key).ToList(); //selecting angles that do not collide
            if ((listY.Count != yAxis.Count) && listY.Count > 0)
            {
                yAngle = CalculateAngle(listY);
                Rotate(cylinder, new Vector3(0, (float)yAngle, 0));
            }

            //outputcylinder();
            //Console.WriteLine("");

            //Third rotation check in Z axis 
            zAxis = CollisionAngles(new Vector3(0, 0, (float)degree)); //creates a dictionary of angles and if the angles causes collision 
            totalAngle = zAxis.Keys.Max();
            Rotate(cylinder, new Vector3(0, 0, -(float)totalAngle)); //rotate back to original position
            List<double> listZ = zAxis.Where(kvp => kvp.Value == false).Select(kvp => kvp.Key).ToList(); //selecting angles that do not collide
            if ((listZ.Count != zAxis.Count) && listZ.Count > 0)
            {
                zAngle = CalculateAngle(listZ);
                Rotate(cylinder, new Vector3(0, 0, (float)zAngle));
            }

            //outputcylinder();
            //Console.WriteLine("");

            //First rotation in X axis
            xAxis = CollisionAngles(new Vector3((float)degree, 0, 0)); //creates a dictionary of angles and if the angles causes collision 
            totalAngle = xAxis.Keys.Max();
            Rotate(cylinder, new Vector3(-(float)totalAngle, 0, 0)); //rotate back to original position
            List<double> listX = xAxis.Where(kvp => kvp.Value == false).Select(kvp => kvp.Key).ToList(); //selecting angles that do not collide
            if ((listX.Count != xAxis.Count) && listX.Count > 0)
            {
                xAngle = CalculateAngle(listX);
                Rotate(cylinder, new Vector3((float)xAngle, 0, 0));
            }

            //Method to find to create a dictionary with angles and if the angles causes collision
            Dictionary<double, bool> CollisionAngles(Vector3 angles)
            {
                Dictionary<double, bool> Axis = new Dictionary<double, bool>();
                for (int i = 1; i <= (180.0 / degree); i++)
                {
                    Rotate(cylinder, angles);
                    totalAngle = i * degree;
                    Axis.Add((totalAngle), CollisionDetect(caList, cylinder));
                }
                return Axis;
            }

            //Method to calculate the angle of rotation
            double CalculateAngle(List<double> AngleList)
            {
                double angleOfRotation = 0;
                if (AngleList[0] == degree && AngleList.Last() == (180 - (180 % degree))) //condition where the axis starts from a non-collision position
                {
                    for (int i = 0; i < AngleList.Count - 1; i++)
                    {

                        if (AngleList[i + 1] - AngleList[i] == degree)
                        {
                            AngleList[i] += 180;
                        }
                        else
                        {
                            AngleList[i] += 180;
                            break;
                        }
                    }
                    //AngleList.ForEach(Console.WriteLine);
                    angleOfRotation = AngleList.Average();
                }
                else
                {
                    angleOfRotation = AngleList.Average();
                }
                return angleOfRotation;
            }
            #endregion

            cylinder.Collision = CollisionDetect(caList, cylinder);
            Console.WriteLine($"Collision is {cylinder.Collision}");
            outputcylinder();
            void outputcylinder()
            {
                cylinder.pointOne += caCentroid;
                cylinder.pointTwo += caCentroid;
                cylinder.centerPoint += caCentroid;
                SharedFunctions.writePymolScriptCylinder(group, cylinder, Global.OUTPUT_DIR, Global.DB_DIR, PdbName);
                cylinder.pointOne -= caCentroid;
                cylinder.pointTwo -= caCentroid;
                cylinder.centerPoint -= caCentroid;
            }

            //Move cylinder back to centroid of the actual protein
            cylinder.pointOne += caCentroid;
            cylinder.pointTwo += caCentroid;
            cylinder.centerPoint += caCentroid;
            cylinder.moment = Vector3.Cross(cylinder.pointOne, cylinder.pointTwo);


            Cylinder Rotate(Cylinder rCylinder, Vector3 rAngles)
            {
                Matrix4x4 rotationMatrixX = new Matrix4x4(1, 0, 0, 0,
                    0, (float)Math.Cos(rAngles.X), (float)Math.Sin(rAngles.X), 0,
                    0, -1 * (float)Math.Sin(rAngles.X), (float)Math.Cos(rAngles.X), 0,
                    0, 0, 0, 1);
                Matrix4x4 rotationMatrixY = new Matrix4x4((float)Math.Cos(rAngles.Y), 0, -1 * (float)Math.Sin(rAngles.Y), 0,
                    0, 1, 0, 0,
                    (float)(Math.Sin(rAngles.Y)), 0, (float)Math.Cos(rAngles.Y), 0,
                    0, 0, 0, 1);
                Matrix4x4 rotationMatrixZ = new Matrix4x4((float)Math.Cos(rAngles.Z), (float)Math.Sin(rAngles.Z), 0, 0,
                    -1 * (float)Math.Sin(rAngles.Z), (float)Math.Cos(rAngles.Z), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);

                rCylinder.pointOne = applyRotation(rCylinder.pointOne, rCylinder.centerPoint);
                rCylinder.pointTwo = applyRotation(rCylinder.pointTwo, rCylinder.centerPoint);

                Vector3 applyRotation(Vector3 point, Vector3 pivot)
                {
                    var pointDiff = point - pivot;
                    pointDiff = Vector3.Transform(pointDiff, rotationMatrixX);
                    pointDiff = Vector3.Transform(pointDiff, rotationMatrixY);
                    pointDiff = Vector3.Transform(pointDiff, rotationMatrixZ);
                    point = pivot + pointDiff;
                    return point;
                }

                return rCylinder;
            }


            bool CollisionDetect(List<Vector3> listOfCA, Cylinder cCylinder)
            {
                var rA = cCylinder.pointOne;
                var rB = cCylinder.pointTwo;
                var e = rA - rB;
                var m = cCylinder.moment;
                var R = cCylinder.radius;

                var collision = listOfCA.AsParallel().Any(ca => ((m + Vector3.Cross(e, ca)).Length()) / e.Length() <= R);

                return collision;
            }

            return cylinder;
        }


        #region Log Data
        private void LogStrandData(GroupOfStrands group)
        {
            string strandLogFileLoc = Global.OUTPUT_DIR + "StrandData.txt";
            using (StreamWriter log = File.AppendText(strandLogFileLoc))
            {
                foreach (var strand in group.StrandSet)
                {
                    if (!strand.Concatenated)
                    {
                        log.Write($"{this.PdbName}\t{strand.StrandInTheGroup.betaStrandNum}\t" + "false" + "\t");
                        foreach (var residue in strand.StrandInTheGroup)
                        {
                            log.Write($"{residue.ThreeLetCode}\t");
                        }
                        log.Write("\n");
                    }
                    else
                    {
                        log.Write($"{this.PdbName}\t{strand.StrandInTheGroup.betaStrandNum}\t" + "true" + "\n");
                    }
                }
            }
        }
        private void LogBarrel(bool IsBarrelFlag, string numberOfStrands)
        {
            string logFileLoc = Global.OUTPUT_DIR + "Log.txt";
            using (StreamWriter log = File.AppendText(logFileLoc))
            {
                log.Write(IsBarrelFlag.ToString() + "\t");
                log.Write(numberOfStrands + "\t");
            }
        }
        #endregion
    }

    public class GroupOfStrands
    {
        public List<OneStrand> StrandSet { get; set; }
        public Coordinate Centroid { get; set; }
        public Cylinder Cylinder { get; set; }
        public bool IsBarrel { get; set; }
        public List<Vector3> NellipseCoords { get; set; }
        public List<Vector3> CellipseCoords { get; set; }
        public Vector3 Ncentroid { get; set; }
        public Vector3 Ccentroid { get; set; }
        public Vector3 Axis { get; set; }
        public GroupOfStrands()
        {
            this.StrandSet = new List<OneStrand>();
            this.Centroid = new Coordinate();
            this.IsBarrel = false;

        }
    }

    public class OneStrand
    {
        public Strand StrandInTheGroup { get; set; }
        public List<Strand> LeftNeighbors { get; set; }
        public List<Strand> RightNeighbors { get; set; }
        public bool DirectionAssigned { get; set; }
        public bool IsAntiParallel { get; set; }
        public bool Concatenated { get; set; }

        public OneStrand(Strand strand)
        {
            this.StrandInTheGroup = strand;

            this.LeftNeighbors = new List<Strand>();
            this.RightNeighbors = new List<Strand>();

            this.DirectionAssigned = false;
            this.IsAntiParallel = false;
            this.Concatenated = false;
        }
    }

    public class Cylinder
    {
        public Vector3 pointOne { get; set; }
        public Vector3 pointTwo { get; set; }

        public Vector3 moment { get; set; }
        public Vector3 centerPoint { get; set; }
        public double radius { get; set; }

        public bool Collision { get; set; }
        public Cylinder()
        {
            this.pointOne = new Vector3(0, 0, 0);
            this.pointTwo = new Vector3(0, 0, 0);
            this.centerPoint = new Vector3(0, 0, 0);
            this.radius = 0;
            this.moment = Vector3.Cross(pointOne, pointTwo);
            this.Collision = true;
        }
    }

}