# Data taken from
# https://en.wikipedia.org/wiki/Electron_configurations_of_the_elements_(data_page)
const table_of_elements = [:H => "hydrogen",
                           :He => "helium",
                           :Li => "lithium",
                           :Be => "beryllium",
                           :B => "boron",
                           :C => "carbon",
                           :N => "nitrogen",
                           :O => "oxygen",
                           :F => "fluorine",
                           :Ne => "neon",
                           :Na => "sodium",
                           :Mg => "magnesium",
                           :Al => "aluminium",
                           :Si => "silicon",
                           :P => "phosphorus",
                           :S => "sulfur",
                           :Cl => "chlorine",
                           :Ar => "argon",
                           :K => "potassium",
                           :Ca => "calcium",
                           :Sc => "scandium",
                           :Ti => "titanium",
                           :V => "vanadium",
                           :Cr => "chromium",
                           :Mn => "manganese",
                           :Fe => "iron",
                           :Co => "cobalt",
                           :Ni => "nickel",
                           :Cu => "copper",
                           :Zn => "zin",
                           :Ga => "gallium",
                           :Ge => "germanium",
                           :As => "arseni",
                           :Se => "selenium",
                           :Br => "bromine",
                           :Kr => "krypton",
                           :Rb => "rubidium",
                           :Sr => "strontium",
                           :Y => "yttrium",
                           :Zr => "zirconium",
                           :Nb => "niobium",
                           :Mo => "molybdenum",
                           :Tc => "technetium",
                           :Ru => "ruthenium",
                           :Rh => "rhodium",
                           :Pd => "palladium",
                           :Ag => "silver",
                           :Cd => "cadmium",
                           :In => "indium",
                           :Sn => "tin",
                           :Sb => "antimony",
                           :Te => "tellurium",
                           :I => "iodine",
                           :Xe => "xenon",
                           :Cs => "caesium",
                           :Ba => "barium",
                           :La => "lanthanum",
                           :Ce => "cerium",
                           :Pr => "praseodymium",
                           :Nd => "neodymium",
                           :Pm => "promethium",
                           :Sm => "samarium",
                           :Eu => "europium",
                           :Gd => "gadolinium",
                           :Tb => "terbium",
                           :Dy => "dysprosium",
                           :Ho => "holmium",
                           :Er => "erbium",
                           :Tm => "thulium",
                           :Yb => "ytterbium",
                           :Lu => "lutetium",
                           :Hf => "hafnium",
                           :Ta => "tantalum",
                           :W => "tungsten",
                           :Re => "rhenium",
                           :Os => "osmium",
                           :Ir => "iridium",
                           :Pt => "platinum",
                           :Au => "gold",
                           :Hg => "mercury",
                           :Tl => "thallium",
                           :Pb => "lead",
                           :Bi => "bismuth",
                           :Po => "polonium",
                           :At => "astatine",
                           :Rn => "radon",
                           :Fr => "francium",
                           :Ra => "radium",
                           :Ac => "actinium",
                           :Th => "thorium",
                           :Pa => "protactinium",
                           :U => "uranium",
                           :Np => "neptunium",
                           :Pu => "plutonium",
                           :Am => "americium",
                           :Cm => "curium",
                           :Bk => "berkelium",
                           :Cf => "californium",
                           :Es => "einsteinium",
                           :Fm => "fermium",
                           :Md => "mendelevium",
                           :No => "nobelium",
                           :Lr => "lawrencium",
                           :Rf => "rutherfordium",
                           :Db => "dubnium",
                           :Sg => "seaborgium",
                           :Bh => "bohrium",
                           :Hs => "hassium",
                           :Mt => "meitnerium",
                           :Ds => "darmstadtium",
                           :Rg => "roentgenium",
                           :Cn => "copernicium",
                           :Nh => "nihonium",
                           :Fl => "flerovium",
                           :Mc => "moscovium",
                           :Lv => "livermorium",
                           :Ts => "tennessine",
                           :Og => "oganesson"]

function element_number(s::Symbol)
    i = findfirst(e -> first(e) == s, table_of_elements)
    i === nothing &&
        throw(ArgumentError("Invalid element $s"))
    i
end
