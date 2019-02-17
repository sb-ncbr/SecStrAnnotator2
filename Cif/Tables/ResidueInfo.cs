namespace Cif.Tables
{
    public struct ResidueInfo {
        public int SeqNumber;
        public string Compound;
        public int? AuthSeqNumber;
        public string AuthInsertionCode;
        public string AuthCompound;

        public ResidueInfo(int seqNumber, string compound, int? authSeqNumber = null, string authInsertionCode = null, string authCompound = null){
            this.SeqNumber = seqNumber;
            this.Compound = compound;
            this.AuthSeqNumber = authSeqNumber;
            this.AuthInsertionCode = authInsertionCode;
            this.AuthCompound = authCompound;
        }
    }
}