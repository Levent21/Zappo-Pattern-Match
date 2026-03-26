# 1. INPUT: Hier kopierst du deine FASTA-Daten einfach zwischen die Anführungszeichen
fasta_input = """
>sp|A0A663DJA2|ORF10_SARS2/1-38 Putative ORF10 protein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=ORF10 PE=5 SV=1
MGYINVFAFPFTIYSLLLCRMNSRNYIAQVDVVNFNLT
>tr|A0A177B238|A0A177B238_9BILA/1-89 Uncharacterized protein (Fragment) OS=Intoshia linei OX=1819745 GN=A3Q56_04062 PE=4 SV=1
MISMYLNLFKFKHLGFHCYNRICHYPIIIYNILLCKFNLRKWYVRVDENVVIGALPINNVLEMLHQTEKVKG
IISLVEKYEMKRFIQFQ
>tr|A0A9K3E1J7|A0A9K3E1J7_HELAN/1-114 Uncharacterized protein OS=Helianthus annuus OX=4232 GN=HanXRQr2_Chr15g0689751 PE=4 SV=1
MKRDALYAVPQKGREAEPMCPATAVKRARASNDKKDRRIRVLKQGGSANGASCFNHIGWENVFNYSISFYRL
GCCNFTGRHADGSLPVEDFNLELENQSGEKVKLGTKDDQRYN
>tr|A0A0F4LB77|A0A0F4LB77_9LACO/1-112 Uncharacterized protein OS=Lactobacillus kimbladii OX=1218506 GN=JF75_16990 PE=4 SV=1
MLKYELNPKIALFLMGILQFLIGILFIYYHKNKIAIYTYFFFGVILLLGAIFIKTTHPVKIIQKANLDGVEK
FCKYLNLALLAFSLLAILFDRTNRTSYIPQIDAFLFFILF
>tr|A0A2H3JIH0|A0A2H3JIH0_WOLCO/1-118 Uncharacterized protein OS=Wolfiporia cocos (strain MD-104) OX=742152 GN=WOLCODRAFT_23438 PE=4 SV=1
MGFHSLLFRMRRMPLLQMFMPLPEGRWLSDASMLVREGELQRAGIIHLMRMGDVVVGDEGNMRRMVWNGGYL
VDLDFSWSCVGDVPQYLPTLAFPPSYFHRIVCTMDSRNPIVHIDIL
"""

# 2. DEFINITION: Die Zappo-Zahlen-Logik (1-7)
ZAPPO_MAP = {
    'I':'1', 'L':'1', 'V':'1', 'A':'1', 'M':'1', # Hautfarbe
    'F':'2', 'W':'2', 'Y':'2',                   # Orange
    'K':'3', 'R':'3', 'H':'3',                   # Blau (Positive)
    'D':'4', 'E':'4',                            # Rot (Negative)
    'S':'5', 'T':'5', 'N':'5', 'Q':'5',          # Grün (Polar)
    'P':'6', 'G':'6',                            # Lila (Spezial)
    'C':'7'                                      # Gelb (Cystein)
}

# 3. OPERATION: Zeile für Zeile durchgehen
def run_zappo_translation(data):
    lines = data.strip().split('\n')
    for line in lines:
        line = line.strip()
        if not line: continue
        
        # Header (mit >) wird 1:1 ausgegeben
        if line.startswith(">"):
            print(line)
        else:
            # Sequenz wird Buchstabe für Buchstabe übersetzt (Zellen-Logik)
            translated = "".join([ZAPPO_MAP.get(aa.upper(), aa) for aa in line])
            print(translated)

# 4. OUTPUT: Starten
run_zappo_translation(fasta_input)

def master_validation(data):
    import re
    print(f"\n{'HEADER-INFO':<20} | {'SOLL (Header)':<12} | {'IST (Zappo)':<10} | {'CHECK'}")
    print("-" * 65)
    
    entries = data.strip().split('>')
    for entry in entries:
        if not entry.strip(): continue
        
        lines = entry.strip().split('\n')
        header = lines[0]
        
        # Sucht nach dem Muster "/1-38" oder "/1-114" im Header
        match = re.search(r'/1-(\d+)', header)
        soll_laenge = int(match.group(1)) if match else "???"
        
        # Die Sequenz zusammenbauen und übersetzen
        full_seq = "".join(lines[1:]).replace(" ", "").replace("\r", "")
        ist_laenge = len(full_seq)
        
        # Vergleich
        status = "✅ PERFECT" if soll_laenge == ist_laenge else "❌ MISMATCH"
        
        # Kurzer Header-Name für die Tabelle
        short_name = header.split('|')[1] if '|' in header else header[:15]
        
        print(f"{short_name:<20} | {soll_laenge:<12} | {ist_laenge:<10} | {status}")

# Aufruf
master_validation(fasta_input)


def analyze_identity_boost(data):
    # 1. Vorbereitung: Wir sammeln die Namen, AS-Sequenzen und Zappo-Ketten
    sequences = []
    entries = data.strip().split('>')
    for entry in entries:
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        as_seq = "".join(lines[1:]).replace(" ", "").upper()
        # Hier nutzen wir deine ZAPPO_MAP für die interne Umwandlung
        zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
        sequences.append({'name': name, 'as': as_seq, 'zap': zap_seq})

    print(f"\n--- MUSTER-VERGLEICH (IDENTITÄTS-CHECK) ---")
    print(f"{'Paar':<30} | {'AS-ID (%)':<10} | {'ZAPPO-ID (%)':<12} | {'BOOST'}")
    print("-" * 70)

    # 2. Vergleich: Jedes Protein mit jedem anderen vergleichen
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            s1, s2 = sequences[i], sequences[j]
            
            # Wir vergleichen auf der Länge des kürzeren Partners
            length = min(len(s1['as']), len(s2['as']))
            
            # Zählen der Treffer
            as_matches = sum(1 for k in range(length) if s1['as'][k] == s2['as'][k])
            zap_matches = sum(1 for k in range(length) if s1['zap'][k] == s2['zap'][k])
            
            # Prozentrechnung
            as_id = (as_matches / length) * 100
            zap_id = (zap_matches / length) * 100
            boost = zap_id - as_id
            
            print(f"{s1['name']} vs {s2['name']:<16} | {as_id:>9.1f} | {zap_id:>12.1f} | +{boost:.1f}%")

# Den neuen Teil jetzt auch aufrufen
analyze_identity_boost(fasta_input)


# --- NEUE SEKTION: SMITH-WATERMAN ALIGNMENT (Lücken-toleranter Vergleich) ---

def run_smith_waterman(seq1, seq2, match=3, mismatch=-1, gap=-2):
    """ Berechnet den Smith-Waterman Score für lokales Alignment """
    n, m = len(seq1), len(seq2)
    # Erstelle die Matrix mit Nullen
    matrix = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Score für Übereinstimmung oder Mismatch
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            # Smith-Waterman Logik (nimmt immer das Maximum, mind. 0)
            matrix[i][j] = max(
                0,
                matrix[i-1][j-1] + score, # Diagonal (Match/Mismatch)
                matrix[i-1][j] + gap,     # Oben (Gap in seq2)
                matrix[i][j-1] + gap      # Links (Gap in seq1)
            )
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
    
    # Normalisierung: Was wäre der perfekte Score der kürzeren Sequenz?
    max_possible = min(n, m) * match
    return (max_score / max_possible) * 100 if max_possible > 0 else 0

def analyze_waterman_boost(data):
    sequences = []
    entries = data.strip().split('>')
    for entry in entries:
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        as_seq = "".join(lines[1:]).replace(" ", "").upper()
        zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
        sequences.append({'name': name, 'as': as_seq, 'zap': zap_seq})

    print(f"\n--- SMITH-WATERMAN VERGLEICH (Lokal mit Gaps) ---")
    print(f"{'Paar':<30} | {'AS-SW %':<10} | {'ZAPPO-SW %':<12} | {'BOOST'}")
    print("-" * 75)

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            s1, s2 = sequences[i], sequences[j]
            
            # 1. SW-Score für Aminosäuren
            as_sw = run_smith_waterman(s1['as'], s2['as'])
            
            # 2. SW-Score für Zappo-Zahlen
            zap_sw = run_smith_waterman(s1['zap'], s2['zap'])
            
            boost = zap_sw - as_sw
            print(f"{s1['name']} vs {s2['name']:<16} | {as_sw:>9.1f} | {zap_sw:>12.1f} | +{boost:.1f}%")

# Start der neuen Analyse
analyze_waterman_boost(fasta_input)


# --- NEUE SEKTION: ZAPPO PATTERN DETECTOR (HSPs) ---

def find_zappo_hsp(data, min_len=5):
    """ Findet die längsten identischen Zahlen-Muster zwischen Sequenzen """
    sequences = []
    entries = data.strip().split('>')
    for entry in entries:
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        as_seq = "".join(lines[1:]).replace(" ", "").upper()
        zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
        sequences.append({'name': name, 'zap': zap_seq})

    print(f"\n--- HIGH-SCORING-PATTERNS (Identische Zappo-Blöcke >= {min_len}) ---")
    print(f"{'Paar':<30} | {'Muster':<15} | {'Länge':<6} | {'Positionen'}")
    print("-" * 75)

    import itertools
    for a, b in itertools.combinations(sequences, 2):
        best_pattern = ""
        pos_a, pos_b = 0, 0
        
        # Suche nach dem längsten gemeinsamen Substring
        for length in range(min_len, min(len(a['zap']), len(b['zap'])) + 1):
            found_in_round = False
            for i in range(len(a['zap']) - length + 1):
                sub = a['zap'][i:i+length]
                if sub in b['zap']:
                    best_pattern = sub
                    pos_a = i + 1
                    pos_b = b['zap'].find(sub) + 1
                    found_in_round = True
            if not found_in_round and best_pattern:
                break
        
        if best_pattern:
            pair_name = f"{a['name']} vs {b['name']}"
            print(f"{pair_name:<30} | {best_pattern:<15} | {len(best_pattern):<6} | A:{pos_a} / B:{pos_b}")

# Aufruf der Muster-Analyse
find_zappo_hsp(fasta_input)


def verify_top_pattern(data, name_a, name_b, pattern):
    sequences = {}
    for entry in data.strip().split('>'):
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        seq = "".join(lines[1:]).replace(" ", "").upper()
        zap = "".join([ZAPPO_MAP.get(aa, '0') for aa in seq])
        sequences[name] = {'as': seq, 'zap': zap}

    s1, s2 = sequences[name_a], sequences[name_b]
    pos1 = s1['zap'].find(pattern)
    pos2 = s2['zap'].find(pattern)

    print(f"\n--- VERIFIZIERUNG: {name_a} vs {name_b} ---")
    print(f"Muster gefunden: {pattern}")
    print(f"{name_a}: {s1['as'][pos1:pos1+len(pattern)]} (Pos {pos1+1})")
    print(f"{name_b}: {s2['as'][pos2:pos2+len(pattern)]} (Pos {pos2+1})")
    print(f"Zappo:  {pattern}")

# Teste es für dein Top-Paar (Pass die Namen evtl. an deine Header an)
verify_top_pattern(fasta_input, "A0A663DJA2", "A0A177B238", "12511173")


# --- ERGÄNZUNG A: FREQUENZ-ANALYSE (2er & 3er Muster) ---
def analyze_pattern_frequency(data):
    sequences = []
    for entry in data.strip().split('>'):
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        as_seq = "".join(lines[1:]).replace(" ", "").upper()
        zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
        sequences.append({'name': name, 'zap': zap_seq, 'as': as_seq})

    print(f"\n--- CHEMISCHER BAUKASTEN (Häufige Muster) ---")
    print(f"{'Typ':<10} | {'Muster':<8} | {'Vorkommen'}")
    print("-" * 40)

    for n in [2, 3, 4, 5, 6, 7, 8]: # Checkt 2er bis 8er Muster
        counts = {}
        for s in sequences:
            seen_in_seq = set([s['zap'][i:i+n] for i in range(len(s['zap'])-(n-1))])
            for pat in seen_in_seq:
                counts[pat] = counts.get(pat, 0) + 1
        
        # Zeige nur Muster, die in fast allen Proteinen (4/5) vorkommen
        for pat, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            if count >= 4:
                label = "Bigramm" if n == 2 else "Trigramm"
                print(f"{label:<10} | {pat:<8} | In {count}/5 Proteinen")

analyze_pattern_frequency(fasta_input)


# --- ERGÄNZUNG: REGIONEN-SCANNER (Hydrophobe vs. Geladene Cluster) ---

def analyze_chemical_clusters(data, window_size=5):
    sequences = []
    entries = data.strip().split('>')
    for entry in entries:
        if not entry.strip(): continue
        lines = entry.strip().split('\n')
        name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
        as_seq = "".join(lines[1:]).replace(" ", "").upper()
        zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
        sequences.append({'name': name, 'zap': zap_seq, 'as': as_seq})

    print(f"\n--- STRUKTURELLE REGIONEN (Cluster-Analyse) ---")
    print(f"{'Protein':<15} | {'Typ':<12} | {'Bereich (AS)':<15} | {'Zappo-Sequenz'}")
    print("-" * 75)

    for s in sequences:
        # Suche nach hydrophoben Clustern (1 und 2)
        for i in range(len(s['zap']) - window_size + 1):
            window = s['zap'][i:i+window_size]
            
            # Check: Besteht das Fenster nur aus 1 oder 2? (Hydrophob/Aromatisch)
            if all(c in '12' for c in window):
                print(f"{s['name']:<15} | HYDROPHOB   | {s['as'][i:i+window_size]:<15} | {window}")
            
            # Check: Besteht das Fenster nur aus 3 oder 4? (Geladen +/-)
            if all(c in '34' for c in window):
                print(f"{s['name']:<15} | GELADEN     | {s['as'][i:i+window_size]:<15} | {window}")

# Aufruf der Cluster-Analyse
analyze_chemical_clusters(fasta_input)


def check_orf10_presence(data, patterns_to_check):
    # ORF10 extrahieren
    entries = data.strip().split('>')
    orf10_zap = ""
    for entry in entries:
        if "ORF10_SARS2" in entry:
            lines = entry.strip().split('\n')
            as_seq = "".join(lines[1:]).replace(" ", "").upper()
            orf10_zap = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
            break
    
    print(f"\n--- ORF10 MUSTER-CHECK ---")
    for pat in patterns_to_check:
        found = "GEFUNDEN ✅" if pat in orf10_zap else "FEHLT ❌"
        pos = orf10_zap.find(pat) + 1 if pat in orf10_zap else "--"
        print(f"Muster {pat:<6} | Status: {found:<12} | Position: {pos}")

# Hier die Muster eintragen, die du in der Liste gesehen hast
check_orf10_presence(fasta_input, ["5111", "1511", "121", "141"])


# --- STARTBEFEHL FÜR DEN DEEP SCAN ---

# 1. Wir müssen die Sequenzen erst vorbereiten (wie in deinen anderen Funktionen)
all_seqs = []
for entry in fasta_input.strip().split('>'):
    if not entry.strip(): continue
    lines = entry.strip().split('\n')
    name = lines[0].split('|')[1] if '|' in lines[0] else lines[0][:15]
    as_seq = "".join(lines[1:]).replace(" ", "").upper()
    zap_seq = "".join([ZAPPO_MAP.get(aa, '0') for aa in as_seq])
    all_seqs.append({'name': name, 'zap': zap_seq})

# 2. JETZT rufen wir die neue Funktion auf
# Wir scannen ORF10 (Index 0) gegen die anderen 4
print("\n" + "="*80)
print("ERGEBNISSE DES DEEP SCANS (2er bis 8er Muster)")
print("="*80)

def orf10_deep_scan_final(sequences, min_len=2, max_len=8):
    target = sequences[0]
    others = sequences[1:]
    all_results = []
    for length in range(min_len, max_len + 1):
        for i in range(len(target['zap']) - length + 1):
            pattern = target['zap'][i:i+length]
            found_in = [s['name'] for s in others if pattern in s['zap']]
            if found_in:
                all_results.append({'pattern': pattern, 'len': length, 'pos': i+1, 'count': len(found_in), 'partners': ", ".join(found_in)})
    
    all_results.sort(key=lambda x: (x['count'], x['len']), reverse=True)
    
    seen = set()
    for r in all_results:
        # Filter: Zeige 2er/3er nur wenn sie in mind. 3 Partnern vorkommen
        # Zeige 4er bis 8er immer, wenn sie mind. 1 Partner haben
        if (r['len'] < 4 and r['count'] >= 3) or (r['len'] >= 4):
            if r['pattern'] not in seen:
                print(f"Muster: {r['pattern']:<10} | L: {r['len']} | Pos: {r['pos']:>2} | Treffer: {r['count']}/4 | Partner: {r['partners']}")
                seen.add(r['pattern'])

# Startet den Scan
orf10_deep_scan_final(all_seqs)
