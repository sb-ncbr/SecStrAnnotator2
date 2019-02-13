namespace Cif.Tables
{
    public struct AtomInfo {
        public string Name;
        public string Element;
        public string AltLoc;
        public bool IsHetatm;
        public double X;
        public double Y;
        public double Z;

        public AtomInfo(string name, string element, string altLoc, bool isHetatm, double x, double y, double z){
            Name = name;
            Element = element;
            AltLoc = altLoc;
            IsHetatm = isHetatm;
            X = x;
            Y = y;
            Z = z;
        }
    }
}