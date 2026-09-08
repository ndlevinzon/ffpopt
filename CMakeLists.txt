cmake_minimum_required(VERSION 3.15)

# =====================================================================
# ACADEMIC ENVIRONMENT CHECK
# =====================================================================
# Industry deployments are constrained from utilizing specific ML model arrays 
# without prior developer coordination. This block checks the parent shell context.
if(DEFINED ENV{ACADEMIC})
    if($ENV{ACADEMIC})
      set(ACADEMIC TRUE)
    else()
      set(ACADEMIC FALSE)
    endif()
  else()
    set(ACADEMIC FALSE)
endif()

if(NOT DEFINED FFPOPTVERSION)
  set(FFPOPTVERSION 0.0)
endif()

# Set model activation variables to TRUE by default
if(NOT DEFINED USE_FENNOL)
  set(USE_FENNOL TRUE)
endif()
if(NOT DEFINED USE_MACE)
  set(USE_MACE TRUE)
endif()
if(NOT DEFINED USE_PM6ML)
  set(USE_PM6ML TRUE)
endif()
if(NOT DEFINED USE_ESPALOMA)
  set(USE_ESPALOMA TRUE)
endif()
if(NOT DEFINED USE_HILFIKER)
  set(USE_HILFIKER TRUE)
endif()

# Configure DFTB parameter sourcing dependencies based on the active AMBER path
if(NOT DEFINED USE_DFTB)
  if(DEFINED ENV{AMBERHOME})
    set(USE_DFTB TRUE)
    message(WARNING "USE_DFTB is TRUE because ENV{AMBERHOME} is defined: $ENV{AMBERHOME}")
  endif()
endif()
if(USE_DFTB)
  if(NOT DEFINED ENV{AMBERHOME})
    set(USE_DFTB FALSE)
    message(WARNING "Disabling DFTB2/DFTB3 parameter downloads because AMBERHOME is unset")
  endif()
endif()

# Apply licensing constraints to target flags if the build is running under non-academic boundaries
if(NOT ACADEMIC)
  set(USE_FENNOL FALSE)
  set(USE_MACE FALSE)
  set(USE_PM6ML FALSE)
endif()

if(NOT USE_MACE)
  set(USE_HILFIKER FALSE)
endif()

project(ffpopt VERSION "${FFPOPTVERSION}" LANGUAGES Fortran)
include(FetchContent)
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development)
set(PackageName "ffpopt") 
set(PKGDATA_DIR "${Python_SITELIB}/${PackageName}/pkgdata")

# Initialize global overrule forcing flags
if(NOT DEFINED FORCE)
  set(FORCE FALSE)
endif()

# FIXED VERIFIED LOCATION: Enforces global HTTP/1.1 serialization down to 
# CMake's cURL engine safely, clearing proxy/REFUSED_STREAM errors cleanly.
set(CMAKE_HTTP_VERSION "1.1")

# Standard base repository identifier to safeguard formatting blocks
set(BASE_URL_DOMAIN "https://github.com")

# =====================================================================
# FFPOPT VERBOSE PERSISTENT REPOSITORY HYPOTHESIS TESTING SUITE
# =====================================================================
message(STATUS "\n" ===========================================================================)
message(STATUS "FFPOPT VERBOSE PERSISTENT REPOSITORY HYPOTHESIS SUITE")
message(STATUS ===========================================================================)
message(STATUS "Forced fresh downloads override (FORCE)        : ${FORCE}")
message(STATUS "Active Shell Terminal Execution Context (PWD)  : $ENV{PWD}")
message(STATUS "---------------------------------------------------------------------------")

# ---------------------------------------------------------------------
# FENNOL MODELS PIPELINE PERCEPTION
# ---------------------------------------------------------------------
if(USE_FENNOL)
  if(FALSE)
    string(JOIN "/" FENNOL_GIT_URL "${BASE_URL_DOMAIN}" "FeNNol-tools" "FeNNol-PMC.git")
    FetchContent_Declare(fennol_pmc GIT_REPOSITORY "${FENNOL_GIT_URL}" GIT_TAG main ENVIRONMENT GIT_LFS_SKIP_SMUDGE=1)
    FetchContent_MakeAvailable(fennol_pmc)
    set(FENNOL_MODEL_DIR "${fennol_pmc_SOURCE_DIR}/FENNIX-BIO1/v1.0")
    set(MODEL_M_PATH "${FENNOL_MODEL_DIR}/fennix-bio1M.fnx")
    set(MODEL_S_PATH "${FENNOL_MODEL_DIR}/fennix-bio1S.fnx")
  else()
    include(ExternalProject)

    # Immutable root path maps pointing directly to your local python code tree folder structure
    set(REPO_M_DIR  "${CMAKE_SOURCE_DIR}/src/python/lib/ffpopt/pkgdata/fennix")
    set(REPO_M_PATH "${REPO_M_DIR}/fennix-bio1M.fnx")
    set(REPO_S_PATH "${REPO_M_DIR}/fennix-bio1S.fnx")

    # Audit file 1: Medium FENNOL parameter set in persistent repository path
    message(STATUS "Probing source tree location for FENNOL Medium binary...")
    message(STATUS "  -> Target Path: ${REPO_M_PATH}")
    if(EXISTS "${REPO_M_PATH}")
      message(STATUS "  -> File Discovered Locally: YES")
      set(EXISTS_M YES)
    else()
      message(STATUS "  -> File Discovered Locally: NO")
      set(EXISTS_M NO)
    endif()

    # Audit file 2: Small FENNOL parameter set in persistent repository path
    message(STATUS "Probing source tree location for FENNOL Small binary...")
    message(STATUS "  -> Target Path: ${REPO_S_PATH}")
    if(EXISTS "${REPO_S_PATH}")
      message(STATUS "  -> File Discovered Locally: YES")
      set(EXISTS_S YES)
    else()
      message(STATUS "  -> File Discovered Locally: NO")
      set(EXISTS_S NO)
    endif()

    # If the files exist inside the source code directory, completely bypass the download step
    if(EXISTS_M AND EXISTS_S AND NOT FORCE)
      set(MODEL_M_PATH "${REPO_M_PATH}")
      set(MODEL_S_PATH "${REPO_S_PATH}")
      set(SKIP_FENNOL_DOWNLOAD TRUE)
      message(STATUS "[DECISION] FENNOL parameters present inside repository tree. Skipping network download.")
    else()
      set(SKIP_FENNOL_DOWNLOAD FALSE)
      string(JOIN "/" FENNOL_M_URL "${BASE_URL_DOMAIN}" "FeNNol-tools" "FeNNol-PMC" "blob" "main" "FENNIX-BIO1" "v1.0" "fennix-bio1M.fnx")
      string(JOIN "/" FENNOL_S_URL "${BASE_URL_DOMAIN}" "FeNNol-tools" "FeNNol-PMC" "blob" "main" "FENNIX-BIO1" "v1.0" "fennix-bio1S.fnx")
      message(STATUS "[DECISION] FENNOL parameters missing from code tree or FORCE active. Initiating network download:")
      message(STATUS "  -> URL Medium: ${FENNOL_M_URL}")
      message(STATUS "  -> URL Small : ${FENNOL_S_URL}")

      # CLEAN EXPLICIT MACROS: Dropped the invalid local HTTP_VERSION parameters to secure positional parsing arrays
      ExternalProject_Add(fennol_model_m
        URL "${FENNOL_M_URL}"
        DOWNLOAD_NO_EXTRACT 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        PREFIX "${CMAKE_BINARY_DIR}/fennixmodels_m"
      )
      ExternalProject_Add(fennol_model_s
        URL "${FENNOL_S_URL}"
        DOWNLOAD_NO_EXTRACT 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        PREFIX "${CMAKE_BINARY_DIR}/fennixmodels_s"
      )
      ExternalProject_Get_Property(fennol_model_m DOWNLOAD_DIR)
      set(MODEL_M_PATH "${DOWNLOAD_DIR}/fennix-bio1M.fnx")
      ExternalProject_Get_Property(fennol_model_s DOWNLOAD_DIR)
      set(MODEL_S_PATH "${DOWNLOAD_DIR}/fennix-bio1S.fnx")
    endif()
  endif()
endif()
# ---------------------------------------------------------------------
# ESPALOMA MODELS PIPELINE PERCEPTION
# ---------------------------------------------------------------------
if(USE_ESPALOMA)
  set(REPO_ESPALOMA_DIR  "${CMAKE_SOURCE_DIR}/src/python/lib/ffpopt/pkgdata/espaloma")
  set(REPO_ESPALOMA_PATH "${REPO_ESPALOMA_DIR}/model.pt")

  message(STATUS "Probing source tree location for ESPALOMA model...")
  message(STATUS "  -> Target Path: ${REPO_ESPALOMA_PATH}")
  if(EXISTS "${REPO_ESPALOMA_PATH}")
    message(STATUS "  -> File Discovered Locally: YES")
    set(EXISTS_ESPALOMA YES)
  else()
    message(STATUS "  -> File Discovered Locally: NO")
    set(EXISTS_ESPALOMA NO)
  endif()

  if(EXISTS_ESPALOMA AND NOT FORCE)
    set(ESPALOMA_PATH "${REPO_ESPALOMA_PATH}")
    set(SKIP_ESPALOMA_DOWNLOAD TRUE)
    message(STATUS "[DECISION] ESPALOMA model present inside repository tree. Skipping network download.")
  else()
    include(ExternalProject)
    set(SKIP_ESPALOMA_DOWNLOAD FALSE)
    string(JOIN "/" ESPALOMA_URL "${BASE_URL_DOMAIN}" "choderalab" "espaloma_charge" "releases" "download" "v0.0.8" "model.pt")
    message(STATUS "[DECISION] ESPALOMA parameters missing or FORCE active. Launching download macro:")
    message(STATUS "  -> Fetching: ${ESPALOMA_URL}")

    ExternalProject_Add(espaloma_model
      URL "${ESPALOMA_URL}"
      DOWNLOAD_NO_EXTRACT 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      PREFIX "${CMAKE_BINARY_DIR}/espaloma_model"
    )
    ExternalProject_Get_Property(espaloma_model DOWNLOAD_DIR)
    set(ESPALOMA_PATH "${DOWNLOAD_DIR}/model.pt")
  endif()
endif()

# ---------------------------------------------------------------------
# HILFIKER MODELS PIPELINE PERCEPTION
# ---------------------------------------------------------------------
if(USE_HILFIKER)
  set(REPO_HILFIKER_DIR  "${CMAKE_SOURCE_DIR}/src/python/lib/ffpopt/pkgdata/ml_for_charges")
  set(REPO_HILFIKER_PATH "${REPO_HILFIKER_DIR}/xgb_chgs.json")

  message(STATUS "Probing source tree location for HILFIKER model...")
  message(STATUS "  -> Target Path: ${REPO_HILFIKER_PATH}")
  if(EXISTS "${REPO_HILFIKER_PATH}")
    message(STATUS "  -> File Discovered Locally: YES")
    set(EXISTS_HILFIKER YES)
  else()
    message(STATUS "  -> File Discovered Locally: NO")
    set(EXISTS_HILFIKER NO)
  endif()

  if(EXISTS_HILFIKER AND NOT FORCE)
    set(HILFIKER_PATH "${REPO_HILFIKER_PATH}")
    set(SKIP_HILFIKER_DOWNLOAD TRUE)
    message(STATUS "[DECISION] HILFIKER model present inside repository tree. Skipping network download.")
  else()
    include(ExternalProject)
    set(SKIP_HILFIKER_DOWNLOAD FALSE)
    string(JOIN "/" HILFIKER_URL "${BASE_URL_DOMAIN}" "mathilfiker" "ml_for_charges" "raw" "refs" "heads" "main" "charge_prediction" "model" "xgb_chgs.json")
    message(STATUS "[DECISION] HILFIKER parameters missing or FORCE active. Launching download macro:")
    message(STATUS "  -> Fetching: ${HILFIKER_URL}")

    ExternalProject_Add(hilfiker_model
      URL "${HILFIKER_URL}"
      DOWNLOAD_NO_EXTRACT 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      PREFIX "${CMAKE_BINARY_DIR}/hilfiker_model"
    )
    ExternalProject_Get_Property(hilfiker_model DOWNLOAD_DIR)
    set(HILFIKER_PATH "${DOWNLOAD_DIR}/xgb_chgs.json")  
  endif()
endif()

message(STATUS ===========================================================================\n)

if(USE_MACE)
  string(JOIN "/" MACE_GIT_URL "${BASE_URL_DOMAIN}" "ACEsuit" "mace-off.git")
  FetchContent_Declare(mace_off GIT_REPOSITORY "${MACE_GIT_URL}" GIT_TAG main)
  FetchContent_MakeAvailable(mace_off)
  set(MACE_MODEL_DIR "${mace_off_SOURCE_DIR}")
endif()

if(USE_PM6ML)
  string(JOIN "/" PM6_GIT_URL "${BASE_URL_DOMAIN}" "Honza-R" "mopac-ml.git")
  FetchContent_Declare(pm6ml GIT_REPOSITORY "${PM6_GIT_URL}" GIT_TAG main)
  FetchContent_MakeAvailable(pm6ml)
  set(PM6ML_MODEL_DIR "${pm6ml_SOURCE_DIR}/models")
endif()

if(USE_DFTB)
  string(JOIN "/" DFTB3_GIT_URL "${BASE_URL_DOMAIN}" "dftbparams" "3ob.git")
  FetchContent_Declare(dftb3obparams GIT_REPOSITORY "${DFTB3_GIT_URL}" GIT_TAG main)
  FetchContent_MakeAvailable(dftb3obparams)
  file(GLOB SRC_FILES "${dftb3obparams_SOURCE_DIR}/skfiles/*.skf")
  set(DFTB3OB_FILES "${SRC_FILES}")

  string(JOIN "/" MIO_GIT_URL "${BASE_URL_DOMAIN}" "dftbparams" "mio.git")
  FetchContent_Declare(dftbmioparams GIT_REPOSITORY "${MIO_GIT_URL}" GIT_TAG main)
  FetchContent_MakeAvailable(dftbmioparams)
  file(GLOB SRC_FILES "${dftbmioparams_SOURCE_DIR}/skfiles/*.skf")
  set(DFTBMIO_FILES "${SRC_FILES}")
endif()

# =====================================================================
# INSTALLATION STEP MANAGEMENT WITH CRITICAL DUAL PACKAGING COPIERS
# =====================================================================
if(DEFINED SKBUILD_SCRIPTS_DIR AND NOT SKBUILD_SCRIPTS_DIR STREQUAL "")
  add_executable(ffpopt-respf src/resp/resp.f)
  set_target_properties(ffpopt-respf PROPERTIES
    Fortran_STANDARD 90
    Fortran_STANDARD_REQUIRED ON
    Fortran_EXTENSIONS OFF
    COMPILE_FLAGS "-std=legacy"
  )
  install(TARGETS ffpopt-respf DESTINATION "${SKBUILD_SCRIPTS_DIR}")
  
  # EXPLANATION: Native install files commands operate outside of pip's build sandbox constraints.
  # If a fresh download happens, it saves copies directly back into your stable git workspace directories.
  if(USE_FENNOL)
    if(NOT SKIP_FENNOL_DOWNLOAD)
      install(FILES "${MODEL_M_PATH}" "${MODEL_S_PATH}" DESTINATION "${REPO_M_DIR}")
    endif()
    install(FILES "${MODEL_M_PATH}" "${MODEL_S_PATH}" DESTINATION "ffpopt/pkgdata/fennix")
  endif()

  if(USE_ESPALOMA)
    if(NOT SKIP_ESPALOMA_DOWNLOAD)
      install(FILES "${ESPALOMA_PATH}" DESTINATION "${REPO_ESPALOMA_DIR}")
    endif()
    install(FILES "${ESPALOMA_PATH}" DESTINATION "ffpopt/pkgdata/espaloma")
  endif()

  if(USE_HILFIKER)
    if(NOT SKIP_HILFIKER_DOWNLOAD)
      install(FILES "${HILFIKER_PATH}" DESTINATION "${REPO_HILFIKER_DIR}")
    endif()
    install(FILES "${HILFIKER_PATH}" DESTINATION "ffpopt/pkgdata/ml_for_charges")
  endif()

  if(USE_MACE)
    install(FILES "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_small.model" "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_medium.model" "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_large.model" DESTINATION "ffpopt/pkgdata/mace-off/mace_off23")
    install(FILES "${MACE_MODEL_DIR}/mace_off24/MACE-OFF24_medium.model" DESTINATION "ffpopt/pkgdata/mace-off/mace_off24")
  endif()
  if(USE_PM6ML)
    install(FILES "${PM6ML_MODEL_DIR}/PM6-ML_correction_seed8_best.ckpt" DESTINATION "ffpopt/pkgdata/pm6ml")
  endif()
  if(USE_DFTB)
    file(COPY ${DFTBMIO_FILES} DESTINATION "$ENV{AMBERHOME}/dat/slko/mio-1-1")
    file(COPY ${DFTB3OB_FILES} DESTINATION "$ENV{AMBERHOME}/dat/slko/3ob-3-1")
  endif()
else()
  add_subdirectory(src/resp)
  if(BUILD_PYTHON)
    add_subdirectory(src/python)
  endif()

  if(USE_FENNOL)
    if(NOT SKIP_FENNOL_DOWNLOAD)
      file(COPY "${MODEL_M_PATH}" "${MODEL_S_PATH}" DESTINATION "${REPO_M_DIR}")
    endif()
    file(COPY "${MODEL_M_PATH}" "${MODEL_S_PATH}" DESTINATION "${PKGDATA_DIR}/fennix")
  endif()

  if(USE_ESPALOMA)
    if(NOT SKIP_ESPALOMA_DOWNLOAD)
      file(COPY "${ESPALOMA_PATH}" DESTINATION "${REPO_ESPALOMA_DIR}")
    endif()
    file(COPY "${ESPALOMA_PATH}" DESTINATION "${PKGDATA_DIR}/espaloma")
  endif()

  if(USE_HILFIKER)
    if(NOT SKIP_HILFIKER_DOWNLOAD)
      file(COPY "${HILFIKER_PATH}" DESTINATION "${REPO_HILFIKER_DIR}")
    endif()
    file(COPY "${HILFIKER_PATH}" DESTINATION "${PKGDATA_DIR}/ml_for_charges")
  endif()

  if(USE_MACE)
    file(COPY "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_small.model" "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_medium.model" "${MACE_MODEL_DIR}/mace_off23/MACE-OFF23_large.model" DESTINATION "${PKGDATA_DIR}/mace-off/mace_off23")
    file(COPY "${MACE_MODEL_DIR}/mace_off24/MACE-OFF24_medium.model" DESTINATION "${PKGDATA_DIR}/mace-off/mace_off24")
  endif()
  if(USE_PM6ML)
    file(COPY "${PM6ML_MODEL_DIR}/PM6-ML_correction_seed8_best.ckpt" DESTINATION "${PKGDATA_DIR}/pm6ml")
  endif()
  if(USE_DFTB)
    file(COPY ${DFTBMIO_FILES} DESTINATION "$ENV{AMBERHOME}/dat/slko/mio-1-1")
    file(COPY ${DFTB3OB_FILES} DESTINATION "$ENV{AMBERHOME}/dat/slko/3ob-3-1")
  endif()
endif()
