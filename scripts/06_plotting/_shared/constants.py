"""Shared constants for plotting topics."""

from __future__ import annotations

TOOL_MAP = {
    "scapa": "scAPA",
    "scapatrap": "scAPAtrap",
    "sierra": "Sierra",
    "maaper": "MAAPER",
    "scapture": "SCAPTURE",
    "scape": "SCAPE",
    "infernape": "Infernape",
    "dapars2": "DaPars2",
    "scmapa": "scMAPA",
}

PROTOCOL_MAP = {
    "Visium": "10X Visium",
    "VisiumHD": "10X Visium HD",
    "Chromium": "10X Chromium",
    "Dropseq": "Drop-seq",
    "Stereoseq": "Stereo-seq",
    "Slideseq": "Slide-seq V2",
    "SpatialTranscriptomics": "ST",
    "Microwell": "Microwell-seq",
    "annotation": "Annotated PAS",
    "Annotation": "Annotated PAS",
    "anno": "Annotated PAS",
}

PROTOCOL_ORDER = [
    "10X Chromium",
    "Drop-seq",
    "Microwell-seq",
    "10X Visium",
    "Stereo-seq",
    "Slide-seq V2",
    "ST",
]

TOOL_ORDER = [
    "scAPAtrap",
    "SCAPE",
    "Infernape",
    "SCAPTURE",
    "scAPA",
    "Sierra",
    "DaPars2",
    "scMAPA",
]

PALETTE_HEX = [
    "#386b98",
    "#269a51",
    "#edaa4d",
    "#d34123",
    "#7e648a",
    "#454545",
    "#929292",
]

TOOL_COLOR_HEX = {
    "scAPAtrap": "#386b98",
    "SCAPE": "#269a51",
    "Infernape": "#edaa4d",
    "SCAPTURE": "#d34123",
    "scAPA": "#7e648a",
    "Sierra": "#454545",
    "DaPars2": "#3b9ab2",
    "scMAPA": "#c76d9e",
}
