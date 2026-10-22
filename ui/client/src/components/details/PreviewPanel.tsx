import { ProjectDocument } from '../../shared/types';
import { useProjectStore } from '../../stores/useProjectStore';
import { ImagePreview, isImageResource } from './ImagePreview';
import { GeoJsonPreview, isGeoJsonResource } from './maps/GeoJsonPreview';
import { isTileUrl, TileMapView } from './maps/TileMapView';

export const hasPreview = (doc: ProjectDocument): boolean => {
  return isTileUrl(doc.resource) || isImageResource(doc.resource) || isGeoJsonResource(doc.resource);
};

export const PreviewPanel = ({
  docid,
}: {
  docid: string;
}) => {
  // The store already holds every document of the project, so there is nothing
  // to fetch: the tab the docid came from is open on a document that is there.
  const project = useProjectStore(state => state.currProject);
  const resource = project?.documents.find(d => d._id.$oid === docid)?.resource;

  if (isTileUrl(resource)) {
    return <TileMapView url={resource} />;
  }
  if (isGeoJsonResource(resource)) {
    return <GeoJsonPreview path={resource} />;
  }
  if (isImageResource(resource)) {
    return <ImagePreview path={resource} />;
  }
  return null;
};
