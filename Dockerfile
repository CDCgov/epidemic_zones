# Multi-platform image (linux/amd64 + linux/arm64).
# Build with:
#   docker buildx build --platform linux/amd64,linux/arm64 -t <name>:<tag> --push .
# Or for a single local platform:
#   docker build -t <name>:<tag> .

FROM rocker/rstudio:4.5.2

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    # HTTP / SSL / XML  (curl, httr, RCurl, xml2, openssl)
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    # Git bindings  (gert)
    libgit2-dev \
    # Spatial  (sf, terra, s2, wk, units)
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libsqlite3-dev \
    libudunits2-dev \
    libabsl-dev \
    cmake \
    # Text / font rendering  (ragg, textshaping, systemfonts, gdtools, svglite, flextable)
    libcairo2-dev \
    libpango1.0-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    # JavaScript engine  (V8)
    libv8-dev \
    # OpenGL / software rendering  (rgl headless via OSMesa)
    libglu1-mesa-dev \
    libgl-dev \
    libosmesa6-dev \
    # ODBC  (odbc)
    unixodbc-dev \
    # Linear programming  (igraph)
    libglpk-dev \
    # Non-linear optimization  (nloptr)
    libnlopt-dev \
    # Tcl/Tk runtime  (tcltk → misc3d, rgl)
    tcl8.6-dev \
    tk8.6-dev \
    # Git required by renv to install GitHub-sourced packages
    git \
    # Pandoc required by rmarkdown / knitr
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Use null OpenGL device so rgl works in headless containers
ENV RGL_USE_NULL=TRUE

# The project .Rprofile sources renv/activate.R, which would repoint the library
# path at the (empty) project-local renv library and hide the packages baked into
# the image. Packages are already pinned to renv.lock versions at build time.
# Written to Renviron.site as well because RStudio Server sessions do not
# reliably inherit container-level ENV.
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
RUN echo 'RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /usr/local/lib/R/etc/Renviron.site \
    && echo 'RGL_USE_NULL=TRUE' >> /usr/local/lib/R/etc/Renviron.site

# ── renv ─────────────────────────────────────────────────────────────────────
# Install the exact renv version recorded in renv.lock (1.1.7)
RUN R -e "install.packages( \
    'https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.7.tar.gz', \
    repos = NULL, type = 'source')"

# ── Project packages ──────────────────────────────────────────────────────────
WORKDIR /project

# Point renv at Posit Package Manager for pre-built Linux binaries.
# renv::restore() still pins every package to the exact version in renv.lock.
ENV RENV_CONFIG_PPM_ENABLED=TRUE

# renv.lock's Repositories field points at an internal CDC-only mirror
# (rpackage.edav.cdc.gov), which is unreachable from outside the CDC network/VPN.
# This override forces renv to fetch from the public Posit Package Manager mirror
# instead, regardless of what's recorded in the lockfile.
ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/latest

COPY renv.lock renv.lock

# Install all 302 packages at the exact versions recorded in renv.lock.
# Packages are written to the default system library so they are available
# in any R session without needing to activate renv.
RUN R -e "renv::restore( \
    lockfile = 'renv.lock', \
    library  = '/usr/local/lib/R/site-library', \
    prompt   = FALSE)" \
    # Strip compiled object/archive files left behind by source builds (~1-3 GB)
    && find /usr/local/lib/R/site-library -name "*.o" -delete \
    && find /usr/local/lib/R/site-library -name "*.a" -delete \
    && rm -rf /tmp/Rtmp*

# Fail the build if any lockfile package is missing (restore reports some
# install failures as warnings only).
RUN R --vanilla -q -e "\
    want <- names(renv::lockfile_read('renv.lock')\$Packages); \
    have <- rownames(installed.packages(lib.loc = .libPaths())); \
    miss <- setdiff(want, have); \
    if (length(miss)) stop('Missing packages: ', paste(miss, collapse = ', ')); \
    cat('All', length(want), 'lockfile packages installed\n')"

# Default credentials: user=rstudio, password set via PASSWORD env var at runtime
ENV PASSWORD=epidemic
EXPOSE 8787
CMD ["/init"]
