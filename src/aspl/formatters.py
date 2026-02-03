from .events import CassetteExon, SpliceJunction, SpliceSite, Exon
from .events import Strand, CoordinateSystem, SpliceSiteType
from dataclasses import dataclass
from typing import TypeVar, cast

T = TypeVar("T", CassetteExon, SpliceJunction, SpliceSite, Exon)


@dataclass
class UnderscoreSeparated:
    coordinate_system: CoordinateSystem = CoordinateSystem.ONE_BASED
    include_site_type: bool = False

    # def __get_number_of_splits(self, target_type: type[T]):
    #     match target_type:
    #         case t if t is CassetteExon:
    #             return (5, 1, )
    #         case t if t is SpliceJunction:
    #             return (3)
    #         case t if t is SpliceSite:
    #             if self.include_site_type:
    #                 return

    def format(self, value: CassetteExon | SpliceJunction | SpliceSite | Exon) -> str:
        assert value.is_valid()
        sites = value.get_sites_sorted_by_coordinate()
        coords = [site.coord for site in sites]
        # add proper checks?
        is_splice_site = len(coords) == 1

        if is_splice_site and self.coordinate_system == CoordinateSystem.ZERO_BASED:
            coords[0] = coords[0] - 1
        fields = [sites[0].seqname] + [str(c) for c in coords] + [sites[0].strand]
        if is_splice_site and self.include_site_type:
            fields += [sites[0].type]
        return "_".join(fields)

    def parse(self, input: str, target_type: type[T]) -> T:
        match target_type:
            case t if t is CassetteExon:
                seqname, *coords, strand = input.rsplit("_", maxsplit=5)
                strand = Strand(strand)
                coords = [int(c) for c in coords]
                if strand == Strand.MINUS:
                    coords = coords[::-1]

                output = CassetteExon(
                    SpliceSite(seqname, coords[0], strand, SpliceSiteType.DONOR),
                    SpliceSite(seqname, coords[1], strand, SpliceSiteType.ACCEPTOR),
                    SpliceSite(seqname, coords[2], strand, SpliceSiteType.DONOR),
                    SpliceSite(seqname, coords[3], strand, SpliceSiteType.ACCEPTOR),
                )

                return cast(T, output)

            case t if t is SpliceJunction:
                seqname, coord1, coord2, strand = input.rsplit("_", maxsplit=3)
                strand = Strand(strand)
                coord1, coord2 = int(coord1), int(coord2)
                if strand == Strand.MINUS:
                    coord2, coord1 = coord1, coord2

                output = SpliceJunction(
                    SpliceSite(seqname, coord1, strand, SpliceSiteType.DONOR),
                    SpliceSite(seqname, coord2, strand, SpliceSiteType.ACCEPTOR),
                )

                return cast(T, output)

            case t if t is SpliceSite:
                if self.include_site_type:
                    seqname, coord, strand, site_type = input.rsplit("_", maxsplit=3)
                    site_type = SpliceSiteType(site_type)
                else:
                    seqname, coord, strand = input.rsplit("_", maxsplit=2)
                    site_type = SpliceSiteType.UNKNOWN

                coord = int(coord)
                if self.coordinate_system == CoordinateSystem.ZERO_BASED:
                    coord += 1

                strand = Strand(strand)
                output = SpliceSite(seqname, coord, strand, site_type)
                return cast(T, output)

            case t if t is CassetteExon:
                raise NotImplementedError("Parsing of exons is not implemented yet")

            case _:
                raise TypeError(f"Unsupported type: {target_type}")
