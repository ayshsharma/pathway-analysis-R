"""
Collect METASPACE datasets that:
1. Use FTICR as analyzer type
2. Have KEGG annotations
3. Belong to selected KEGG species

Install:
    pip install metaspace2020 pandas tqdm

Optional, for private datasets:
    from metaspace import SMInstance
    sm = SMInstance()
    sm.save_login()
"""

from metaspace import SMInstance
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import re



# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

SPECIES = [
    "Homo sapiens (human)",
    "Homo sapiens (human) ",
    "Pan troglodytes",
    "Macaca mulatta",
    "Macaca fascicularis",
    "Mus musculus (mouse)",
    "Rattus norvegicus (rat)",
    "Oryctolagus cuniculus",
    "Canis lupus familiaris",
    "Felis catus",
    "Bos taurus",
    "Bos indicus",
    "Capra hircus",
    "Ovis aries",
    "Sus scrofa",
    "Equus caballus",
    "Equus asinus",
    "Hipposideros armiger",
    "Desmodus rotundus",
    "Canis familiaris",
    "Bovine",
    "Mouse",
    "mouse",
    "Monkey",
    "Rhesus (Monkey)",
    "Non-human primate (rhesus macaque)",
    "Dog",
    "pig skin",
    "Pig",
    "Sus scrofa domesticus (pig)",
    "Gallus domesticus (chicken)",
    "Mouse ",
    "Mouse Brain",
    "Cat",
    "cat",
    "pig (KVL)",
    "Drosophila melanogaster",
    "Apis mellifera",
    "Galleria mallonella",
    "Glycine max (soybean)",
    "Schistosoma mansoni",
    "Maize",
    "Chlamydomonas reinhardtii",
    "Arabidopsis thaliana",
    "Gossypium hirsutum",
]


ANALYZER_TYPE = "FTICR"
FDR_THRESHOLD = 0.50

OUTPUT_DATASETS_CSV = "metaspace_fticr_kegg_datasets.csv"
OUTPUT_ANNOTATIONS_CSV = "metaspace_fticr_kegg_annotations_summary.csv"


ANNOTATION_CSV_DIR = Path("metaspace_kegg_annotation_csvs")
ANNOTATION_CSV_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def safe_getattr(obj, attr, default=None):
    """Safely get an attribute from an object."""
    return getattr(obj, attr, default)


def get_dataset_metadata_value(ds, key, default=None):
    """
    Safely retrieve metadata values if ds.metadata is available.
    METASPACE metadata structures can vary between datasets.
    """
    metadata = safe_getattr(ds, "metadata", None)

    if metadata is None:
        return default

    if isinstance(metadata, dict):
        return metadata.get(key, default)

    return default


def db_value(db, key, default=None):
    """
    Read database fields from either:
    - dict-like objects: db["name"] / db.get("name")
    - object-like entries: db.name, db.version, db.id
    """
    if isinstance(db, dict):
        return db.get(key, default)

    return getattr(db, key, default)


def find_kegg_databases(ds):
    """
    Return KEGG-like database entries from ds.database_details.

    Handles both dict and MolecularDB object entries.
    """
    kegg_dbs = []

    try:
        database_details = ds.database_details
    except Exception as exc:
        print(f"Could not read database details for {ds}: {exc}")
        return kegg_dbs

    for db in database_details:
        db_name = str(db_value(db, "name", ""))
        db_version = str(db_value(db, "version", ""))
        db_id = str(db_value(db, "id", ""))

        searchable = f"{db_name} {db_version} {db_id}".upper()

        if "KEGG" in searchable:
            kegg_dbs.append(db)

    return kegg_dbs


def count_kegg_annotations(ds, kegg_db, fdr_threshold=0.10):
    """
    Count KEGG annotations for a dataset/database.

    Handles both dict and MolecularDB object database entries.
    """
    db_name = db_value(kegg_db, "name")
    db_version = db_value(kegg_db, "version")

    if not db_name:
        return None, None

    try:
        if db_version:
            results = ds.results(database=(db_name, db_version))
        else:
            results = ds.results(database=db_name)

        if results is None or len(results) == 0:
            return 0, 0

        total_annotations = len(results)

        if "fdr" in results.columns:
            annotations_below_fdr = int((results["fdr"] <= fdr_threshold).sum())
        else:
            annotations_below_fdr = None

        return total_annotations, annotations_below_fdr

    except Exception as exc:
        print(
            f"Could not fetch results for dataset {safe_getattr(ds, 'id', ds)} "
            f"and database {db_name} {db_version}: {exc}"
        )
        return None, None


# ---------------------------------------------------------------------
# Main collection logic
# ---------------------------------------------------------------------

def collect_metaspace_fticr_kegg_datasets():
    sm = SMInstance()

    dataset_rows = []
    annotation_rows = []

    seen_dataset_ids = set()

    for species in SPECIES:
        print(f"\nSearching species: {species}")

        try:
            datasets = sm.datasets(
                analyzer_type=ANALYZER_TYPE,
                organism=species,
            )
        except Exception as exc:
            print(f"Search failed for species {species}: {exc}")
            continue

        print(f"Found {len(datasets)} FTICR datasets before KEGG filtering")

        for ds in tqdm(datasets, desc=f"Filtering {species}"):

            ds_id = safe_getattr(ds, "id", None)
            ds_name = safe_getattr(ds, "name", None)

            # Avoid duplicate datasets if the API returns overlap
            unique_key = ds_id or ds_name
            if unique_key in seen_dataset_ids:
                continue

            kegg_dbs = find_kegg_databases(ds)

            if not kegg_dbs:
                continue

            seen_dataset_ids.add(unique_key)

            dataset_row = {
                "dataset_id": ds_id,
                "dataset_name": ds_name,
                "organism_query": species,
                "analyzer_type_query": ANALYZER_TYPE,
                "polarity": safe_getattr(ds, "polarity", None),
                "ionisation_source": safe_getattr(ds, "ionisation_source", None),
                "submitter": safe_getattr(ds, "submitter", None),
                "group": safe_getattr(ds, "group", None),
                "project": safe_getattr(ds, "project", None),
                "url": (
                    f"https://metaspace2020.org/annotations?ds={ds_id}"
                    if ds_id else None
                ),
                "kegg_database_count": len(kegg_dbs),
                
                "kegg_databases": "; ".join(
                    [
                        f"{db_value(db, 'name', '')} {db_value(db, 'version', '')}".strip()
                        for db in kegg_dbs
                    ]
                ),
            }

            dataset_rows.append(dataset_row)

            for kegg_db in kegg_dbs:
                annotation_csv_path, total_annotations, annotations_below_fdr = save_kegg_annotations_csv(
                    ds=ds,
                    kegg_db=kegg_db,
                    organism_query=species,
                    output_dir=ANNOTATION_CSV_DIR,
                )

                annotation_rows.append(
                    {
                        "dataset_id": ds_id,
                        "dataset_name": ds_name,
                        "organism_query": species,
                        "database_id": db_value(kegg_db, "id"),
                        "database_name": db_value(kegg_db, "name"),
                        "database_version": db_value(kegg_db, "version"),
                        "total_annotations": total_annotations,
                        f"annotations_fdr_lte_{FDR_THRESHOLD}": annotations_below_fdr,
                        "annotation_csv_path": annotation_csv_path,
                    }
                )

    datasets_df = pd.DataFrame(dataset_rows)
    annotations_df = pd.DataFrame(annotation_rows)

    datasets_df.to_csv(OUTPUT_DATASETS_CSV, index=False)
    annotations_df.to_csv(OUTPUT_ANNOTATIONS_CSV, index=False)

    print("\nDone.")
    print(f"Datasets written to: {OUTPUT_DATASETS_CSV}")
    print(f"Annotation summary written to: {OUTPUT_ANNOTATIONS_CSV}")
    print(f"Matched datasets: {len(datasets_df)}")

    return datasets_df, annotations_df


def sanitize_filename(value, max_length=120):
    """
    Make a safe filename from dataset/database strings.
    """
    value = str(value)
    value = re.sub(r"[^\w\-.]+", "_", value)
    value = value.strip("_")

    if len(value) > max_length:
        value = value[:max_length]

    return value or "unknown"


def get_kegg_results_dataframe(ds, kegg_db):
    """
    Fetch the actual annotation results table for one dataset and one KEGG database.

    Returns:
        pandas.DataFrame or None
    """
    db_name = db_value(kegg_db, "name")
    db_version = db_value(kegg_db, "version")

    if not db_name:
        return None

    try:
        if db_version:
            results = ds.results(database=(db_name, db_version))
        else:
            results = ds.results(database=db_name)

        if results is None or len(results) == 0:
            return pd.DataFrame()

        # ds.results() often returns formula/adduct as index.
        # reset_index() preserves them as normal CSV columns.
        results_df = results.reset_index()

        return results_df

    except Exception as exc:
        print(
            f"Could not fetch annotation results for dataset "
            f"{safe_getattr(ds, 'id', ds)} and database {db_name} {db_version}: {exc}"
        )
        return None


def save_kegg_annotations_csv(ds, kegg_db, organism_query, output_dir=ANNOTATION_CSV_DIR):
    """
    Save actual KEGG annotation results for one dataset/database to a CSV file.

    Returns:
        annotation_csv_path, total_annotations, annotations_below_fdr
    """
    ds_id = safe_getattr(ds, "id", None)
    ds_name = safe_getattr(ds, "name", None)

    db_id = db_value(kegg_db, "id")
    db_name = db_value(kegg_db, "name")
    db_version = db_value(kegg_db, "version")

    results_df = get_kegg_results_dataframe(ds, kegg_db)

    if results_df is None:
        return None, None, None

    total_annotations = len(results_df)

    if total_annotations == 0:
        annotations_below_fdr = 0
    elif "fdr" in results_df.columns:
        annotations_below_fdr = int((results_df["fdr"] <= FDR_THRESHOLD).sum())
    else:
        annotations_below_fdr = None

    # Add helpful metadata columns to every exported annotation CSV
    results_df.insert(0, "dataset_id", ds_id)
    results_df.insert(1, "dataset_name", ds_name)
    results_df.insert(2, "organism_query", organism_query)
    results_df.insert(3, "database_id", db_id)
    results_df.insert(4, "database_name", db_name)
    results_df.insert(5, "database_version", db_version)

    filename = (
        f"{sanitize_filename(ds_id)}__"
        f"{sanitize_filename(ds_name)}__"
        f"{sanitize_filename(db_name)}__"
        f"{sanitize_filename(db_version)}.csv"
    )

    output_path = output_dir / filename

    results_df.to_csv(output_path, index=False)

    return str(output_path), total_annotations, annotations_below_fdr


if __name__ == "__main__":
    datasets_df, annotations_df = collect_metaspace_fticr_kegg_datasets()

    print("\nDataset preview:")
    print(datasets_df.head())

    print("\nAnnotation summary preview:")
    print(annotations_df.head())