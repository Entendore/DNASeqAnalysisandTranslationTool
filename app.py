import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import zipfile
from io import BytesIO
import io
import pandas as pd
from base64 import b64encode
from fpdf import FPDF

st.set_page_config(layout="wide")
st.title("DNA/RNA Sequence Analysis & Six-Frame ORF Finder with Viewer, Filters & Exports")

# --------- Helper Functions ---------

def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return 100 * gc / len(seq) if len(seq) > 0 else 0

def find_orfs_with_coords(dna_seq: Seq, min_len: int, table: int = 1):
    seq_len = len(dna_seq)
    results = []
    rc_seq = dna_seq.reverse_complement()

    frames = {
        "+1": dna_seq[0:],
        "+2": dna_seq[1:],
        "+3": dna_seq[2:],
        "-1": rc_seq[0:],
        "-2": rc_seq[1:],
        "-3": rc_seq[2:]
    }

    for frame, seq_fragment in frames.items():
        prot_seq = seq_fragment.translate(table=table, to_stop=False)
        aa_start = None
        for i, aa in enumerate(prot_seq):
            if aa != "*" and aa_start is None:
                aa_start = i
            elif aa == "*" and aa_start is not None:
                orf_len = i - aa_start
                if orf_len >= min_len:
                    if frame.startswith("+"):
                        frame_offset = int(frame[1]) - 1
                        dna_start = frame_offset + aa_start * 3
                        dna_end = frame_offset + i * 3 - 1
                        strand = "+"
                        dna_seq_str = str(dna_seq[dna_start:dna_end + 1])
                    else:
                        frame_offset = int(frame[1]) - 1
                        dna_start = seq_len - (frame_offset + i * 3)
                        dna_end = seq_len - (frame_offset + aa_start * 3) - 1
                        strand = "-"
                        dna_seq_str = str(dna_seq[dna_start:dna_end + 1].reverse_complement())
                    prot_str = str(prot_seq[aa_start:i])
                    results.append({
                        "frame": frame,
                        "strand": strand,
                        "prot": prot_str,
                        "dna": dna_seq_str,
                        "start": dna_start + 1,
                        "end": dna_end + 1
                    })
                aa_start = None
        # Handle ORF running to end (no stop codon)
        if aa_start is not None and (len(prot_seq) - aa_start) >= min_len:
            i = len(prot_seq)
            if frame.startswith("+"):
                frame_offset = int(frame[1]) - 1
                dna_start = frame_offset + aa_start * 3
                dna_end = frame_offset + i * 3 - 1
                strand = "+"
                dna_seq_str = str(dna_seq[dna_start:dna_end + 1])
            else:
                frame_offset = int(frame[1]) - 1
                dna_start = seq_len - (frame_offset + i * 3)
                dna_end = seq_len - (frame_offset + aa_start * 3) - 1
                strand = "-"
                dna_seq_str = str(dna_seq[dna_start:dna_end + 1].reverse_complement())
            prot_str = str(prot_seq[aa_start:i])
            results.append({
                "frame": frame,
                "strand": strand,
                "prot": prot_str,
                "dna": dna_seq_str,
                "start": dna_start + 1,
                "end": dna_end + 1
            })

    return results

def make_fai(fasta_str):
    fai_lines = []
    offset = 0
    seq_name = None
    seq_len = 0
    line_bases = None
    line_width = None

    lines = fasta_str.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(">"):
            if seq_name is not None:
                fai_lines.append(f"{seq_name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_width}")
            seq_name = line[1:].split()[0]
            seq_offset = offset
            seq_len = 0
            line_bases = None
            line_width = None
        else:
            seq_len += len(line)
            if line_bases is None:
                line_bases = len(line)
                line_width = len(line) + 1
        offset += len(line) + 1

    if seq_name is not None:
        fai_lines.append(f"{seq_name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_width}")

    return "\n".join(fai_lines) + "\n"

def export_gff3(records_orfs, source="seqtool"):
    gff_lines = ["##gff-version 3"]
    for rec_id, orfs in records_orfs:
        for idx, orf in enumerate(orfs, start=1):
            orf_id = f"{rec_id}_frame{orf['frame']}_orf{idx}_pos{orf['start']}-{orf['end']}_strand{orf['strand']}"
            start = orf["start"]
            end = orf["end"]
            strand = orf["strand"]
            frame_abs = abs(int(orf["frame"][1]))
            phase = (frame_abs - 1) % 3
            attr_str = (
                f"ID={orf_id};"
                f"Name={orf_id};"
                f"protein_id={orf_id};"
                f"Note=Length {len(orf['prot'])} aa"
            )
            gff_lines.append(
                f"{rec_id}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attr_str}"
            )
    return "\n".join(gff_lines) + "\n"

def orf_color(length_aa):
    if length_aa > 150:
        return "green"
    elif length_aa > 50:
        return "orange"
    else:
        return "red"

def highlight_orfs_with_anchors(seq, orfs):
    seq_list = list(seq)
    anchor_positions = {orf['start'] - 1: idx + 1 for idx, orf in enumerate(orfs)}
    highlight_mask = [False] * len(seq_list)

    for orf in orfs:
        start, end = orf['start'] - 1, orf['end']
        for i in range(start, end):
            if 0 <= i < len(seq_list):
                highlight_mask[i] = True

    result = []
    for i, base in enumerate(seq_list):
        if i in anchor_positions:
            result.append(f'<a id="orf{anchor_positions[i]}"></a>')
        if highlight_mask[i]:
            result.append(f"<mark>{base}</mark>")
        else:
            result.append(base)
    return "".join(result)

def get_highlighted_html(seq_id, highlighted_seq_html):
    html = f"""
    <html>
    <head>
    <style>
    body {{ font-family: monospace; white-space: pre-wrap; line-height: 1.2; }}
    mark {{ background-color: yellow; }}
    </style>
    </head>
    <body>
    <h1>Sequence: {seq_id}</h1>
    <pre>{highlighted_seq_html}</pre>
    </body>
    </html>
    """
    return html

def generate_pdf(seq_id, filtered_orfs, seq_str):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Courier", size=10)
    pdf.cell(0, 10, f"Annotated Sequence: {seq_id}", ln=True, align="L")
    pdf.ln(5)

    line_len = 80
    seq_len = len(seq_str)
    orf_positions = [(orf['start'] - 1, orf['end']) for orf in filtered_orfs]

    for i in range(0, seq_len, line_len):
        line_seq = seq_str[i:i+line_len]
        line_display = ""
        for j, base in enumerate(line_seq):
            pos = i + j
            in_orf = any(start <= pos < end for start, end in orf_positions)
            if in_orf:
                line_display += f"[{base.upper()}]"
            else:
                line_display += base.lower()
        pdf.multi_cell(0, 6, line_display)

    pdf.ln(10)
    pdf.cell(0, 10, "ORF Details:", ln=True)
    for idx, orf in enumerate(filtered_orfs, start=1):
        pdf.cell(0, 8,
            f"ORF {idx}: Frame {orf['frame']}, Strand {orf['strand']}, Positions {orf['start']}-{orf['end']}, Length {len(orf['prot'])} aa",
            ln=True)

    pdf_output = BytesIO()
    pdf.output(pdf_output)
    return pdf_output.getvalue()

# --------- MAIN APP ---------

uploaded_file = st.file_uploader("Upload FASTA or FASTQ file", type=["fasta", "fa", "fastq", "fq"])
orf_min_len = st.number_input("Minimum ORF length (aa)", min_value=10, max_value=1000, value=30, step=1)
translate_table = st.number_input("Translation Table (Genetic Code)", min_value=1, max_value=33, value=1, step=1)

if uploaded_file:
    text_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    file_type = "fasta" if uploaded_file.name.lower().endswith(("fa", "fasta")) else "fastq"
    records = list(SeqIO.parse(text_io, file_type))

    if not records:
        st.error("No sequences found.")
    else:
        st.success(f"Loaded {len(records)} sequences.")

        protein_fasta_entries = []
        dna_fasta_entries = []
        orf_rows = []
        codon_counts_total = Counter()
        aa_counts_total = Counter()
        records_orfs = []

        # Find ORFs and gather stats
        for rec in records:
            seq_str = str(rec.seq).upper()
            dna_seq = Seq(seq_str)
            orfs = find_orfs_with_coords(dna_seq, min_len=orf_min_len, table=translate_table)
            records_orfs.append((rec.id, orfs))

            for idx, orf in enumerate(orfs, start=1):
                orf_id = f"{rec.id}_frame{orf['frame']}_orf{idx}_pos{orf['start']}-{orf['end']}_strand{orf['strand']}"
                protein_fasta_entries.append(f">{orf_id}\n{orf['prot']}")
                dna_fasta_entries.append(f">{orf_id}\n{orf['dna']}")
                orf_rows.append({
                    "SeqID": rec.id,
                    "ORF_ID": orf_id,
                    "Frame": orf["frame"],
                    "Strand": orf["strand"],
                    "Start": orf["start"],
                    "End": orf["end"],
                    "Length_nt": len(orf["dna"]),
                    "Length_aa": len(orf["prot"]),
                    "Protein_Seq": orf["prot"]
                })
                for i in range(0, len(orf["dna"]) - 2, 3):
                    codon = orf["dna"][i:i + 3].upper()
                    if len(codon) == 3:
                        codon_counts_total[codon] += 1
                aa_counts_total.update(orf["prot"])

        # Summary table per record
        summary_df = pd.DataFrame([{
            "ID": rec.id,
            "Description": rec.description,
            "Length": len(str(rec.seq)),
            "ORF_Count": len(orfs)
        } for rec, (_, orfs) in zip(records, records_orfs)])

        st.dataframe(summary_df)

        # Filtering UI
        min_aa_length = st.slider("Filter ORFs by minimum length (aa)", 10, 300, orf_min_len)
        strand_filter = st.multiselect("Select strands", options=["+", "-"], default=["+", "-"])
        frame_filter = st.multiselect("Select frames",
                                     options=["+1", "+2", "+3", "-1", "-2", "-3"],
                                     default=["+1", "+2", "+3", "-1", "-2", "-3"])

        seq_choice = st.selectbox("Choose sequence to view", [rec.id for rec in records])
        chosen_seq = next(rec for rec in records if rec.id == seq_choice)
        chosen_orfs = next(orfs for rid, orfs in records_orfs if rid == seq_choice)

        filtered_orfs = [
            orf for orf in chosen_orfs
            if len(orf['prot']) >= min_aa_length
            and orf['strand'] in strand_filter
            and orf['frame'] in frame_filter
        ]

        highlighted_seq = highlight_orfs_with_anchors(str(chosen_seq.seq), filtered_orfs)

        st.markdown(f"### Sequence: {seq_choice} (length {len(chosen_seq.seq)})")
        st.markdown(
            "<div style='font-family: monospace; white-space: pre-wrap; line-height: 1.2; max-height: 400px; overflow-y: scroll;'>"
            f"{highlighted_seq}</div>",
            unsafe_allow_html=True,
        )

        st.markdown(f"### Filtered ORFs ({len(filtered_orfs)})")
        for idx, orf in enumerate(filtered_orfs, start=1):
            color = orf_color(len(orf['prot']))
            st.markdown(
                f'<a href="#orf{idx}" style="color:{color}; font-weight:bold;">ORF {idx}:</a> '
                f'Frame {orf["frame"]}, Strand {orf["strand"]}, Positions {orf["start"]}-{orf["end"]} ({len(orf["prot"])} aa)',
                unsafe_allow_html=True
            )
            st.code(orf['prot'])

        # Inject JS for smooth scrolling
        smooth_scroll_js = """
        <script>
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
        </script>
        """
        components.html(smooth_scroll_js, height=0, width=0)

        # Export files as ZIP
        protein_fasta_str = "\n".join(protein_fasta_entries)
        dna_fasta_str = "\n".join(dna_fasta_entries)
        orf_df = pd.DataFrame(orf_rows)
        orf_csv = orf_df.to_csv(index=False)
        orf_tsv = orf_df.to_csv(index=False, sep="\t")

        codon_df = pd.DataFrame(sorted(codon_counts_total.items()), columns=["Codon", "Count"])
        aa_df = pd.DataFrame(sorted(aa_counts_total.items()), columns=["AminoAcid", "Count"])
        codon_csv = codon_df.to_csv(index=False)
        aa_csv = aa_df.to_csv(index=False)

        protein_fai = make_fai(protein_fasta_str)
        dna_fai = make_fai(dna_fasta_str)
        gff_str = export_gff3(records_orfs)

        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("proteins.faa", protein_fasta_str)
            zf.writestr("proteins.faa.fai", protein_fai)
            zf.writestr("orfs_dna.fna", dna_fasta_str)
            zf.writestr("orfs_dna.fna.fai", dna_fai)
            zf.writestr("orfs.gff3", gff_str)
            zf.writestr("orfs_table.csv", orf_csv)
            zf.writestr("orfs_table.tsv", orf_tsv)
            zf.writestr("codon_usage.csv", codon_csv)
            zf.writestr("aa_composition.csv", aa_csv)

        st.download_button(
            "Download All Results (ZIP)",
            data=zip_buffer.getvalue(),
            file_name="annotation_results.zip",
            mime="application/zip"
        )

        # Export annotated sequence as HTML
        html_str = get_highlighted_html(seq_choice, highlighted_seq)
        b64 = b64encode(html_str.encode()).decode()
        href = f'<a href="data:text/html;base64,{b64}" download="{seq_choice}_annotated.html">Download annotated sequence (HTML)</a>'
        st.markdown(href, unsafe_allow_html=True)

        # PDF export button
        if st.button("Download annotated sequence (PDF)"):
            pdf_bytes = generate_pdf(seq_choice, filtered_orfs, str(chosen_seq.seq))
            st.download_button(
                "Click here to download PDF",
                data=pdf_bytes,
                file_name=f"{seq_choice}_annotated.pdf",
                mime="application/pdf"
            )
