#import "paste.h"
#import <Foundation/Foundation.h>
#import <Cocoa/Cocoa.h>
#import <unistd.h>
#include <algorithm>

bool paste(unsigned char *& IM,size_t & h, size_t & w, size_t & c)
{
  NSPasteboard *pasteboard = [NSPasteboard generalPasteboard];
  NSArray *classArray = [NSArray arrayWithObject:[NSImage class]];
  NSDictionary *options = [NSDictionary dictionary];
  BOOL ok = [pasteboard canReadObjectForClasses:classArray options:options]; 
  if(!ok)
  {
    //printf("Error: clipboard doesn't seem to contain image...");
    return false;
  }
  NSArray *objectsToPaste = [pasteboard readObjectsForClasses:classArray options:options];
  NSImage *image = [objectsToPaste objectAtIndex:0];
  NSBitmapImageRep* bitmap = [NSBitmapImageRep imageRepWithData:[image TIFFRepresentation]];
  // http://stackoverflow.com/a/19649616/148668
  w = [bitmap pixelsWide];
  h = [bitmap pixelsHigh];
  size_t rowBytes = [bitmap bytesPerRow];
  c = rowBytes/w;
  unsigned char* pixels = [bitmap bitmapData];

  //IM = new unsigned char[h*rowBytes];
  //std::copy(pixels,pixels+h*rowBytes,IM);
  //std::swap(h,c);
  //IM = permute(IM,[3 2 1]);

  IM = new unsigned char[w*h*c];
  for(size_t y = 0; y < h ; y++)
  {
    for(size_t x = 0; x < w; x++)
    {
      for(size_t d = 0;d<c;d++)
      {
        // For some reason colors are invertex
        IM[y+h*(x+d*w)] = pixels[y*rowBytes + x*c + d];
      }
    }
  }

  //[image release];
  return true;
}

#ifdef MAIN
int main(int argc, char * const argv[])
{
  unsigned char * IM;
  size_t h,w,c;
  return paste(IM,h,w,c)?EXIT_SUCCESS:EXIT_FAILURE;
}
#endif
