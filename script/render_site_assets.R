# Regenerate every PNG the workflowr site displays *from* the panel PDFs in
# figures/, so that the site can never drift from the panel set the figure
# scripts produced.
#
# Run from the repository root after any script/Figure*.R re-run and before
# workflowr::wflow_build():
#
#     Rscript script/render_site_assets.R
#
# Why convert rather than re-plot: if the site PNGs had their own ggsave()
# calls, any panel whose layout is not fully determined by the data --
# ggrepel labels, jitter, force-directed layouts, the random subsampling
# several panels use -- would render differently in the PNG than in the PDF,
# and the two would diverge silently as scripts were re-run at different
# times. Converting the PDF guarantees the page shows the panel itself.
#
# Resolution: 144 dpi, i.e. 2 px per PDF point, which stays sharp on
# high-density displays. Override with DENSITY=72 for smaller files.
#
# Conversion is done by Ghostscript rather than by an R package, because it is
# the one rasteriser available on every machine this has been run on: magick's
# image_read_pdf() needs pdftools, which needs poppler-cpp headers to build,
# and magick's own PDF delegate needs Ghostscript anyway.

gs_bin <- Sys.which("gs")
if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') not found on PATH; install it (brew install ",
         "ghostscript) or convert figures/**/*.pdf to ",
         "analysis/assets/<Figure>/<panel>.png by other means.")
}

density <- as.numeric(Sys.getenv("DENSITY", "144"))
src <- "figures"
dst <- "analysis/assets"

# "Figure 4" -> "Figure4", "Extended Data Figure 11" -> "ExtendedDataFigure11".
asset_dir_name <- function(x) gsub(" ", "", x)

pdfs <- list.files(src, pattern = "\\.pdf$", recursive = TRUE,
                   full.names = TRUE)
if (!length(pdfs)) stop("No panel PDFs found under ", src, "/")

n <- 0L
for (pdf in pdfs) {
    figure <- basename(dirname(pdf))
    panel <- tools::file_path_sans_ext(basename(pdf))
    outdir <- file.path(dst, asset_dir_name(figure))
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(outdir, paste0(panel, ".png"))

    status <- system2(gs_bin, c(
        "-q", "-dNOPAUSE", "-dBATCH", "-dSAFER",
        "-sDEVICE=png16m", sprintf("-r%g", density),
        "-dFirstPage=1", "-dLastPage=1",
        "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
        shQuote(sprintf("-sOutputFile=%s", out)),
        shQuote(pdf)
    ), stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(out)) {
        stop("Ghostscript failed to convert ", pdf)
    }

    cat(sprintf("  %-58s -> %s\n", pdf, out))
    n <- n + 1L
}
cat(sprintf("converted %d panel PDFs at %g dpi\n", n, density))

# workflowr renders each page with rmarkdown::render() in its own session,
# which does not perform rmarkdown's site-resource copying. docs/assets is
# therefore never refreshed by wflow_build(), not even with republish = TRUE.
# Mirror it here, or the built site keeps serving whatever PNGs were there
# before.
if (dir.exists("docs")) {
    dir.create("docs/assets", recursive = TRUE, showWarnings = FALSE)
    unlink(list.files("docs/assets", full.names = TRUE), recursive = TRUE)
    file.copy(list.files(dst, full.names = TRUE), "docs/assets",
              recursive = TRUE)
    cat(sprintf("mirrored %s -> docs/assets\n", dst))
}

# Fail loudly if a page references a PNG that was not produced.
rmds <- list.files("analysis", pattern = "\\.Rmd$", full.names = TRUE)
refs <- unique(unlist(lapply(rmds, function(f) {
    txt <- readLines(f, warn = FALSE)
    regmatches(txt, gregexpr('assets/[^"]+\\.png', txt))
})))
missing <- refs[!file.exists(file.path("analysis", refs))]
if (length(missing)) {
    cat("MISSING assets referenced by analysis/*.Rmd:\n")
    cat(paste0("  ", missing, collapse = "\n"), "\n")
    stop("site references ", length(missing), " PNG(s) that do not exist")
}
cat("all", length(refs), "include_graphics() references resolve\n")
