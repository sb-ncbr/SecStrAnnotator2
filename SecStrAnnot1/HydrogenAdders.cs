﻿using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;
using Cif.Tables;

namespace protein
{
	public static class HydrogenAdders
	{
		public interface IHydrogenAdder {
			// List<Residue> AddHydrogens (IEnumerable<Residue> residues);
			Protein AddHydrogens (Protein protein);
		}

		public class DsspLikeAmideHydrogenAdder : IHydrogenAdder {
			const double NH_BOND_LENGTH = 1.0;

			public DsspLikeAmideHydrogenAdder(){
			}

			/**
			 * First removes all hydrogens and then adds amide hydrogens using simple approach described in DSSP.
			 * (Position of H is 1.0 Angstrom from N, in the opposite direction than carboxylic O is from C.)
			 */
			public Protein AddHydrogens(Protein protein){
				//TODO implement this somehow!
				ModelBuilder builder = new ModelBuilder();
				foreach (Entity entity in protein.GetEntities()){
					builder.StartEntity(entity.Id);
					foreach (Chain chain in entity.GetChains()) {
						builder.StartChain(chain.Id, chain.AuthId);
						foreach (Fragment fragment in chain.GetFragments()) {
							Residue[] residues = fragment.GetResidues().ToArray();
							for (int i = 0; i < residues.Length; i++) {
								Residue residue = residues[i];
								builder.StartResidue(residue.SeqNumber, residue.Compound);
								foreach (Atom atom in residue.GetAtoms()) {
									builder.AddAtom(atom.Id, atom.AtomInfo());
								}
								if (i > 0){
									AtomInfo? newHydrogen = CalculateNewHydrogen(residues[i-1], residue);
									if (newHydrogen != null){
										builder.AddAtom((AtomInfo) newHydrogen);
									}
								}
							}
						}
					}
				}
				Protein result = new Protein(builder.GetModel(protein.Model.ModelNumber));
				Cif.Libraries.Lib.WriteLineDebug($"AddHydrogens(): {result.Model.Residues.Count}");
				return result;
				// return new Protein(builder.GetModel(protein.Model.ModelNumber));
			}

			private static AtomInfo? CalculateNewHydrogen(Residue previousResidue, Residue thisResidue){
				if (thisResidue.IsProline()){
					return null;
				}
				Atom? nAtom = thisResidue.GetNAmide();
				Atom? cAtom = previousResidue.GetCCarb();
				Atom? oAtom = previousResidue.GetOCarb();
				if (nAtom==null || cAtom==null || oAtom==null){
					return null;
				}
				Vector n = ((Atom)nAtom).Position();
				Vector c = ((Atom)cAtom).Position();
				Vector o = ((Atom)oAtom).Position();
				Vector h = n + NH_BOND_LENGTH * (c - o).Normalize ();

				return new AtomInfo(Atom.NAME_H_AMIDE, Atom.ELEMENT_H, ((Atom)nAtom).AltLoc, ((Atom)nAtom).IsHetatm, h.X, h.Y, h.Z );
			}
		
		}



	}
}

