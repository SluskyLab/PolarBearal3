using System;
using System.Collections.Generic;
using System.Numerics;



namespace betaBarrelProgram
{
    public class Atom
    {
        public string AtomType { get; set; }
        public string AtomName { get; set; }
        public Vector3 Coords { get; set; }
        public int AtomNum { get; set; }
        public int ResSeqID { get; set; }
        public Vector3 Hydrogen { get; set; }
        public Vector3 e1 { get; set; }
        public Vector3 e2 { get; set; }
        public double partialcharge { get; set; }
        public List<Atom> SCSCNeighAtoms { get; set; }
        public List<Atom> SCBBNeighAtoms { get; set; }
        public List<Atom> BBNeighAtoms { get; set; }

        public Atom(string resName, Vector3 coords, int atomNum, string atomName, string atomType)
        {
            this.AtomNum = atomNum;
            this.Coords = coords;
            this.AtomName = atomName;
            this.AtomType = atomType;

            Tuple<string, string> key = new Tuple<string, string>(resName, atomName);
            if (Global.partialChargesDict.ContainsKey(key) == true) this.partialcharge = Global.partialChargesDict[key];
            else if (atomName == "H") this.partialcharge = 0.31;
            else if (atomName == "N") this.partialcharge = -0.47;
            else if (atomName == "C") this.partialcharge = 0.51;
            else if (atomName == "O") this.partialcharge = -0.51;
            else this.partialcharge = 0;

            this.SCSCNeighAtoms = new List<Atom>();
            this.SCBBNeighAtoms = new List<Atom>();
            this.BBNeighAtoms = new List<Atom>();

        }

        public void translate(Vector3 translationVec)
        {

            this.Coords += translationVec;

        }
    }
}